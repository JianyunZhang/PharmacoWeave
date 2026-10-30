"""Single-dose IV bolus NCA.

The curve uses quantified concentrations only. BLQ points stay censored.
AUC is linear-up/log-down. For a profile that does not already contain a
time-zero concentration, C0 is the log-linear back-extrapolation through the
first two quantified points, and the interval from 0 to the first sample is
included in AUC.

The terminal slope uses the adjusted R² rule
``pharmacoweave.lambda_z.adjusted_r2.v1``:

- candidate windows are the last k quantified points at or after Tmax
- k runs from the configured minimum (at least 3) through all such points
- each window is an ordinary least-squares fit of ln(concentration) on time
- the window with the largest adjusted R² is kept, and ties prefer the longer window
- the phase is qualified only when the slope is negative and adjusted R²
  meets the configured minimum

Half-life, AUC∞, CL, and Vz are reported only when that phase is qualified.
An extrapolated AUC fraction above the configured maximum is stored and
warned about. ``auc_inf_is_reportable`` is the reporting flag for that cutoff.
"""

import math
from dataclasses import dataclass

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import NcaParameterSet, ParameterSemantic
from pharmacoweave.core.profiles import ConcentrationProfile
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_nca.config import NcaConfig

TERMINAL_RULE = "pharmacoweave.lambda_z.adjusted_r2.v1"


@dataclass(frozen=True)
class _TerminalPhase:
    qualified: bool
    n_points: int
    adjusted_r2: float | None
    lambda_per_h: float | None


def nca_iv_bolus(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Calculate NCA for one IV bolus profile."""

    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("nca_iv_bolus requires route iv_bolus.")
    warnings = [f"Terminal phase rule: {TERMINAL_RULE}."]
    if profile.n_predose_omitted:
        warnings.append(f"Left {profile.n_predose_omitted} predose rows out of the post-dose curve.")
    if profile.n_blq_omitted or profile.n_above_uloq_omitted:
        warnings.append(
            f"Omitted {profile.n_blq_omitted} BLQ and {profile.n_above_uloq_omitted} above-ULOQ rows. "
            "Those rows were not replaced with a number."
        )
    if not profile.points:
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

    cmax_point = max(profile.points, key=lambda point: (point.concentration_mg_l, -point.time_h))
    curve, c0_warning = _curve_with_optional_c0(profile.points)
    if c0_warning:
        warnings.append(c0_warning)
    auc_last = _auc_linear_up_log_down(curve)
    terminal_points = [point for point in profile.points if point.time_h + 1e-12 >= cmax_point.time_h]
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
        clearance_semantic=None if clearance is None else ParameterSemantic.CL,
        clearance=clearance,
        clearance_unit=None if clearance is None else "L/h",
        volume_semantic=None if volume is None else ParameterSemantic.VZ,
        volume=volume,
        volume_unit=None if volume is None else "L",
        warnings=warnings,
    )


def _curve_with_optional_c0(points: tuple) -> tuple[list[tuple[float, float]], str | None]:
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


def _terminal_phase(points: list, *, min_points: int, min_adjusted_r2: float) -> _TerminalPhase:
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


def _log_regression(points: list) -> tuple[float, float, float] | None:
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
