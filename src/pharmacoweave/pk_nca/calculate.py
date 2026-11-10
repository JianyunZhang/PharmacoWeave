"""Single-dose rich-serial NCA for IV bolus, constant-rate IV infusion, and extravascular routes.

The curve uses quantified concentrations only. BLQ points stay censored.
AUC is linear-up/log-down:

- ``iv_bolus``: when no time-zero sample is present and the first two quantified
  points decline, C0 is back-extrapolated log-linearly and included in AUC;
  candidate terminal windows use points with ``t >= Tmax``; clearance and
  terminal volume are reported as ``CL`` and ``Vz``.
- ``oral``, ``subcutaneous``, ``intramuscular``: integrates from ``(0, 0)``
  without C0 back-extrapolation; candidate terminal windows use strictly
  post-peak points (``t > Tmax``); apparent clearance and terminal volume are
  reported as ``CL_over_F`` and ``Vz_over_F`` because bioavailability ``F`` is
  unknown in a single-arm extravascular design.
- ``iv_infusion``: requires positive ``infusion_duration_h``; integrates from
  ``(0, 0)`` without bolus C0 back-extrapolation; candidate terminal windows
  use post-infusion points with ``t >= max(Tmax, infusion_duration_h)``;
  clearance and terminal volume are reported as ``CL`` and ``Vz``.

The terminal slope uses the adjusted R² rule
``pharmacoweave.lambda_z.adjusted_r2.v1``.
"""

import math
from dataclasses import dataclass

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import NcaParameterSet, ParameterSemantic
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_nca.config import NcaConfig

TERMINAL_RULE = "pharmacoweave.lambda_z.adjusted_r2.v1"
_EXTRAVASCULAR_ROUTES = {Route.ORAL, Route.SUBCUTANEOUS, Route.INTRAMUSCULAR}


@dataclass(frozen=True)
class _TerminalPhase:
    qualified: bool
    n_points: int
    adjusted_r2: float | None
    lambda_per_h: float | None


def nca_profile(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Dispatch NCA according to the config dosing and ``profile.route``."""

    if config.dosing in {"repeated_equal", "steady_state"}:
        return nca_interval(profile, config)
    if profile.route is Route.IV_BOLUS:
        return nca_iv_bolus(profile, config)
    if profile.route is Route.IV_INFUSION:
        return nca_iv_infusion(profile, config)
    if profile.route in _EXTRAVASCULAR_ROUTES:
        return nca_extravascular(profile, config)
    raise AnalysisError(f"NCA is not implemented for route {profile.route.value}.")


def nca_interval(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """NCA over one equal-dose interdose interval.

    AUC_tau is the linear-up/log-down integral only when quantified samples
    exist at both ends. No concentration is imputed. Clearance is Dose/AUC_tau
    only when the config declares ``steady_state``.
    """

    if config.dosing not in {"repeated_equal", "steady_state"}:
        raise AnalysisError("nca_interval requires repeated_equal or steady_state.")
    if profile.tau_h is None or profile.tau_h <= 0:
        raise AnalysisError("Interval NCA requires a positive dosing interval.")
    tau = profile.tau_h
    warnings = _initial_warnings(profile)
    warnings.append(
        "Equal-dose, equal-interval NCA uses observed samples at the interval ends. "
        "Missing ends are not imputed."
    )
    if config.dosing == "repeated_equal":
        warnings.append("Clearance is omitted because Dose/AUC_tau equals clearance only at steady state.")
    else:
        warnings.append(
            "Steady state is declared by the caller. Clearance is Dose/AUC_tau under that declaration, "
            "not a demonstration that steady state was reached."
        )
    if not profile.points:
        return _empty_interval_result(profile, config, warnings, tau)

    snapped = _snap_interval_points(profile.points, tau)
    cmax_point = max(snapped, key=lambda point: (point.concentration_mg_l, -point.time_h))
    curve = [(point.time_h, point.concentration_mg_l) for point in snapped]
    auc_observed = _auc_linear_up_log_down(curve)
    has_start = any(point.time_h == 0.0 for point in snapped)
    has_end = any(point.time_h == tau for point in snapped)
    auc_tau = auc_observed if has_start and has_end else None
    if auc_tau is None:
        warnings.append("AUC_tau is omitted because a quantified sample is missing at 0 h or at tau.")
    cmin = next((point.concentration_mg_l for point in snapped if point.time_h == tau), None)
    cavg = None if auc_tau is None else auc_tau / tau
    fluctuation = None
    if cavg is not None and cmin is not None and cavg > 0:
        fluctuation = (cmax_point.concentration_mg_l - cmin) / cavg
    terminal_points = _interval_terminal_points(profile, snapped, cmax_point)
    terminal = _terminal_phase(
        terminal_points,
        min_points=config.lambda_z_min_points,
        min_adjusted_r2=config.lambda_z_min_adjusted_r2,
    )
    half_life = None
    lambda_z = None
    if terminal.qualified and terminal.lambda_per_h is not None:
        lambda_z = terminal.lambda_per_h
        half_life = math.log(2.0) / lambda_z
    elif terminal.n_points < config.lambda_z_min_points:
        warnings.append("The terminal phase has fewer quantified points than lambda_z_min_points.")
    else:
        warnings.append("The terminal phase did not meet the adjusted R² minimum, so λz is omitted.")
    warnings.append("AUC∞ is not reported for an interdose interval.")
    clearance_semantic = (
        ParameterSemantic.CL_OVER_F if profile.route in _EXTRAVASCULAR_ROUTES else ParameterSemantic.CL
    )
    volume_semantic = (
        ParameterSemantic.VZ_OVER_F if profile.route in _EXTRAVASCULAR_ROUTES else ParameterSemantic.VZ
    )
    clearance = None
    volume = None
    if config.dosing == "steady_state" and auc_tau is not None and auc_tau > 0 and profile.dose_mg > 0:
        clearance = profile.dose_mg / auc_tau
        if lambda_z is not None and lambda_z > 0:
            volume = clearance / lambda_z
    accumulation = None
    if (
        profile.reference_cmax_mg_l is not None
        and profile.reference_cmax_mg_l > 0
        and config.dosing == "repeated_equal"
    ):
        accumulation = cmax_point.concentration_mg_l / profile.reference_cmax_mg_l
    return NcaParameterSet(
        subject_id=profile.subject_id,
        analyte_id=profile.analyte_id,
        route=profile.route,
        dose_event_id=profile.dose_event_id,
        auc_method=config.auc_method,
        terminal_phase_qualified=terminal.qualified,
        auc_last=auc_observed,
        auc_last_unit="mg*h/L",
        cmax=cmax_point.concentration_mg_l,
        cmax_unit="mg/L",
        tmax=cmax_point.time_h,
        tmax_unit="h",
        lambda_z=lambda_z,
        lambda_z_unit=None if lambda_z is None else "1/h",
        half_life=half_life,
        half_life_unit=None if half_life is None else "h",
        terminal_phase_n_points=terminal.n_points,
        terminal_phase_adjusted_r2=terminal.adjusted_r2,
        clearance_semantic=None if clearance is None else clearance_semantic,
        clearance=clearance,
        clearance_unit=None if clearance is None else "L/h",
        volume_semantic=None if volume is None else volume_semantic,
        volume=volume,
        volume_unit=None if volume is None else "L",
        regimen=config.dosing,
        n_doses=profile.n_doses,
        dosing_interval_h=tau,
        auc_tau=auc_tau,
        auc_tau_unit=None if auc_tau is None else "mg*h/L",
        cmin=cmin,
        cmin_unit=None if cmin is None else "mg/L",
        cavg=cavg,
        cavg_unit=None if cavg is None else "mg/L",
        fluctuation=fluctuation,
        accumulation_ratio=accumulation,
        warnings=warnings,
    )


def _empty_interval_result(
    profile: ConcentrationProfile,
    config: NcaConfig,
    warnings: list[str],
    tau: float,
) -> NcaParameterSet:
    warnings.append("No quantified concentrations inside the dosing interval.")
    return NcaParameterSet(
        subject_id=profile.subject_id,
        analyte_id=profile.analyte_id,
        route=profile.route,
        dose_event_id=profile.dose_event_id,
        auc_method=config.auc_method,
        terminal_phase_qualified=False,
        excluded_reason="No quantified concentrations inside the dosing interval.",
        regimen=config.dosing,
        n_doses=profile.n_doses,
        dosing_interval_h=tau,
        warnings=warnings,
    )


def _snap_interval_points(points: tuple[QuantifiedPoint, ...], tau: float) -> list[QuantifiedPoint]:
    snapped: list[QuantifiedPoint] = []
    for point in points:
        if point.time_h < -1e-6 or point.time_h > tau + 1e-6:
            continue
        time_h = point.time_h
        if abs(time_h) <= 1e-6:
            time_h = 0.0
        elif abs(time_h - tau) <= 1e-6:
            time_h = tau
        snapped.append(
            QuantifiedPoint(
                observation_id=point.observation_id,
                time_h=time_h,
                concentration_mg_l=point.concentration_mg_l,
            )
        )
    return snapped


def _interval_terminal_points(
    profile: ConcentrationProfile,
    points: list[QuantifiedPoint],
    cmax_point: QuantifiedPoint,
) -> list[QuantifiedPoint]:
    if profile.route is Route.IV_BOLUS:
        return [point for point in points if point.time_h + 1e-12 >= cmax_point.time_h]
    if profile.route is Route.IV_INFUSION:
        cutoff = cmax_point.time_h
        if profile.infusion_duration_h is not None:
            cutoff = max(cutoff, profile.infusion_duration_h)
        return [point for point in points if point.time_h + 1e-12 >= cutoff]
    return [point for point in points if point.time_h > cmax_point.time_h + 1e-12]


def nca_iv_bolus(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Calculate NCA for one IV bolus profile."""

    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("nca_iv_bolus requires route iv_bolus.")
    warnings = _initial_warnings(profile)
    if not profile.points:
        return _empty_nca_result(profile, config, warnings)

    cmax_point = max(profile.points, key=lambda point: (point.concentration_mg_l, -point.time_h))
    curve, c0_warning = _curve_with_optional_c0(profile.points)
    if c0_warning:
        warnings.append(c0_warning)
    terminal_points = [point for point in profile.points if point.time_h + 1e-12 >= cmax_point.time_h]
    return _finalize_nca(
        profile=profile,
        config=config,
        curve=curve,
        cmax_point=cmax_point,
        terminal_points=terminal_points,
        clearance_semantic=ParameterSemantic.CL,
        volume_semantic=ParameterSemantic.VZ,
        warnings=warnings,
    )


def nca_extravascular(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Calculate NCA for one single-dose extravascular profile (oral, SC, or IM)."""

    if profile.route not in _EXTRAVASCULAR_ROUTES:
        raise AnalysisError("nca_extravascular requires route oral, subcutaneous, or intramuscular.")
    warnings = _initial_warnings(profile)
    if not profile.points:
        return _empty_nca_result(profile, config, warnings)

    cmax_point = max(profile.points, key=lambda point: (point.concentration_mg_l, -point.time_h))
    curve, zero_note = _curve_from_zero(profile.points)
    if zero_note:
        warnings.append(zero_note)
    terminal_points = [point for point in profile.points if point.time_h > cmax_point.time_h + 1e-12]
    return _finalize_nca(
        profile=profile,
        config=config,
        curve=curve,
        cmax_point=cmax_point,
        terminal_points=terminal_points,
        clearance_semantic=ParameterSemantic.CL_OVER_F,
        volume_semantic=ParameterSemantic.VZ_OVER_F,
        warnings=warnings,
    )


def nca_iv_infusion(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Calculate NCA for one single-dose constant-rate IV infusion profile."""

    if profile.route is not Route.IV_INFUSION:
        raise AnalysisError("nca_iv_infusion requires route iv_infusion.")
    if profile.infusion_duration_h is None or profile.infusion_duration_h <= 0:
        raise AnalysisError("nca_iv_infusion requires a positive infusion_duration.")
    warnings = _initial_warnings(profile)
    if not profile.points:
        return _empty_nca_result(profile, config, warnings)

    cmax_point = max(profile.points, key=lambda point: (point.concentration_mg_l, -point.time_h))
    curve, zero_note = _curve_from_zero(profile.points)
    if zero_note:
        warnings.append(zero_note)
    if profile.points[0].time_h > profile.infusion_duration_h + 1e-12:
        warnings.append(
            "First quantified sample occurs after infusion end; "
            "linear trapezoid from (0, 0) may understate infusion AUC."
        )
    cutoff_h = max(cmax_point.time_h, profile.infusion_duration_h)
    terminal_points = [point for point in profile.points if point.time_h + 1e-12 >= cutoff_h]
    return _finalize_nca(
        profile=profile,
        config=config,
        curve=curve,
        cmax_point=cmax_point,
        terminal_points=terminal_points,
        clearance_semantic=ParameterSemantic.CL,
        volume_semantic=ParameterSemantic.VZ,
        warnings=warnings,
    )


def _initial_warnings(profile: ConcentrationProfile) -> list[str]:
    warnings = [f"Terminal phase rule: {TERMINAL_RULE}."]
    if profile.n_predose_omitted:
        warnings.append(f"Left {profile.n_predose_omitted} predose rows out of the post-dose curve.")
    if profile.n_blq_omitted or profile.n_above_uloq_omitted:
        warnings.append(
            f"Omitted {profile.n_blq_omitted} BLQ and {profile.n_above_uloq_omitted} above-ULOQ rows. "
            "Those rows were not replaced with a number."
        )
    return warnings


def _empty_nca_result(
    profile: ConcentrationProfile,
    config: NcaConfig,
    warnings: list[str],
) -> NcaParameterSet:
    return NcaParameterSet(
        subject_id=profile.subject_id,
        analyte_id=profile.analyte_id,
        route=profile.route,
        dose_event_id=profile.dose_event_id,
        auc_method=config.auc_method,
        terminal_phase_qualified=False,
        excluded_reason="No quantified post-dose concentrations.",
        warnings=warnings,
    )


def _finalize_nca(
    *,
    profile: ConcentrationProfile,
    config: NcaConfig,
    curve: list[tuple[float, float]],
    cmax_point: QuantifiedPoint,
    terminal_points: list[QuantifiedPoint],
    clearance_semantic: ParameterSemantic,
    volume_semantic: ParameterSemantic,
    warnings: list[str],
) -> NcaParameterSet:
    auc_last = _auc_linear_up_log_down(curve)
    terminal = _terminal_phase(
        terminal_points,
        min_points=config.lambda_z_min_points,
        min_adjusted_r2=config.lambda_z_min_adjusted_r2,
    )
    auc_inf = None
    fraction = None
    half_life = None
    clearance = None
    volume = None
    if terminal.qualified and terminal.lambda_per_h is not None:
        clast = curve[-1][1]
        extrapolated = clast / terminal.lambda_per_h
        auc_inf = auc_last + extrapolated
        fraction = extrapolated / auc_inf if auc_inf > 0 else None
        half_life = math.log(2.0) / terminal.lambda_per_h
        if profile.dose_mg > 0 and auc_inf > 0:
            clearance = profile.dose_mg / auc_inf
            volume = clearance / terminal.lambda_per_h
        if fraction is not None and fraction > config.auc_inf_max_extrapolated_fraction:
            warnings.append(
                f"AUC∞ extrapolated fraction {fraction:.3f} exceeds "
                f"{config.auc_inf_max_extrapolated_fraction:.2f}. "
                "auc_inf_is_reportable is false for this cutoff."
            )
    elif terminal.n_points < config.lambda_z_min_points:
        warnings.append("The terminal phase has fewer quantified points than lambda_z_min_points.")
    else:
        warnings.append("The terminal phase did not meet the adjusted R² minimum, so λz and AUC∞ are omitted.")

    return NcaParameterSet(
        subject_id=profile.subject_id,
        analyte_id=profile.analyte_id,
        route=profile.route,
        dose_event_id=profile.dose_event_id,
        auc_method=config.auc_method,
        terminal_phase_qualified=terminal.qualified,
        auc_last=auc_last,
        auc_last_unit="mg*h/L",
        auc_inf=auc_inf,
        auc_inf_unit=None if auc_inf is None else "mg*h/L",
        auc_inf_extrap_fraction=fraction,
        cmax=cmax_point.concentration_mg_l,
        cmax_unit="mg/L",
        tmax=cmax_point.time_h,
        tmax_unit="h",
        lambda_z=terminal.lambda_per_h if terminal.qualified else None,
        lambda_z_unit="1/h" if terminal.qualified else None,
        half_life=half_life,
        half_life_unit="h" if half_life is not None else None,
        terminal_phase_n_points=terminal.n_points,
        terminal_phase_adjusted_r2=terminal.adjusted_r2,
        clearance_semantic=None if clearance is None else clearance_semantic,
        clearance=clearance,
        clearance_unit=None if clearance is None else "L/h",
        volume_semantic=None if volume is None else volume_semantic,
        volume=volume,
        volume_unit=None if volume is None else "L",
        warnings=warnings,
    )


def _curve_with_optional_c0(points: tuple[QuantifiedPoint, ...]) -> tuple[list[tuple[float, float]], str | None]:
    series = [(point.time_h, point.concentration_mg_l) for point in points]
    if series[0][0] <= 1e-12:
        return series, None
    if len(series) < 2 or not (series[0][1] > series[1][1] > 0):
        return series, "AUC starts at the first quantified sample because C0 was not back-extrapolated."
    t1, c1 = series[0]
    t2, c2 = series[1]
    slope = math.log(c2 / c1) / (t2 - t1)
    c0 = math.exp(math.log(c1) - slope * t1)
    return [(0.0, c0), *series], (
        "C0 was back-extrapolated on the log-linear line through the first two quantified points and included in AUC."
    )


def _curve_from_zero(points: tuple[QuantifiedPoint, ...]) -> tuple[list[tuple[float, float]], str | None]:
    series = [(point.time_h, point.concentration_mg_l) for point in points]
    if series[0][0] <= 1e-12:
        return series, None
    return [(0.0, 0.0), *series], (
        "AUC integrates from (0, 0) to the first quantified sample without bolus C0 back-extrapolation."
    )


def _auc_linear_up_log_down(series: list[tuple[float, float]]) -> float:
    total = 0.0
    for (t1, c1), (t2, c2) in zip(series, series[1:], strict=False):
        dt = t2 - t1
        if dt < 0:
            raise AnalysisError("Concentration times are not sorted.")
        if dt == 0:
            continue
        if c1 <= 0 or c2 <= 0 or c2 >= c1:
            total += dt * (c1 + c2) / 2.0
        else:
            total += dt * (c2 - c1) / math.log(c2 / c1)
    return total


def _terminal_phase(
    points: list[QuantifiedPoint],
    *,
    min_points: int,
    min_adjusted_r2: float,
) -> _TerminalPhase:
    usable = [point for point in points if point.concentration_mg_l > 0]
    if len(usable) < min_points:
        return _TerminalPhase(False, len(usable), None, None)
    best: tuple[float, int, float] | None = None
    for k in range(min_points, len(usable) + 1):
        window = usable[-k:]
        intercept_slope = _log_regression(window)
        if intercept_slope is None:
            continue
        _intercept, slope, adjusted = intercept_slope
        if slope >= 0:
            continue
        if best is None or adjusted > best[0] or (adjusted == best[0] and k > best[1]):
            best = (adjusted, k, -slope)
    if best is None:
        return _TerminalPhase(False, len(usable), None, None)
    adjusted, k, lambda_per_h = best
    qualified = adjusted >= min_adjusted_r2
    return _TerminalPhase(qualified, k, adjusted, lambda_per_h if qualified else None)


def _log_regression(points: list[QuantifiedPoint]) -> tuple[float, float, float] | None:
    xs = [point.time_h for point in points]
    ys = [math.log(point.concentration_mg_l) for point in points]
    n = len(xs)
    if n < 3:
        return None
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0:
        return None
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    intercept = ybar - slope * xbar
    sst = sum((y - ybar) ** 2 for y in ys)
    sse = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys, strict=True))
    sse = max(0.0, sse)
    if sst <= 1e-18:
        adjusted = 1.0 if sse <= 1e-18 else 0.0
    else:
        adjusted = 1.0 - (sse / sst) * (n - 1) / (n - 2)
    return intercept, slope, adjusted
