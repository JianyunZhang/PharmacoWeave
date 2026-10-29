"""Single-dose IV bolus NCA for rich serial sampling.

The trapezoid is linear on a rising or flat interval and logarithmic on a
declining interval where both concentrations are positive. For an IV bolus
whose first sample is after time zero, C0 is back-extrapolated from the log
slope of the first two declining quantified points and is not reported as Cmax.

The terminal slope uses quantified points strictly after Tmax. Regressions of
the last k points, k from the configured minimum upward, are ordinary least
squares of ln(concentration) on time. The regression with the largest adjusted
R² is kept. Adjusted R² values within 1e-6 are tied, and the tie goes to the
regression with more points. AUC∞, clearance, and Vz are omitted when the
extrapolated fraction exceeds the configured reporting limit. λz and t½ are
omitted when the terminal regression itself is not qualified.

BLQ samples stay censored. They are not given a number and they are not used
as zero. This module does not call the one-compartment closed form.
"""

import math
from dataclasses import dataclass

from pharmacoweave.core.entities import ObservationStatus, Route
from pharmacoweave.core.modeling import NcaParameterSet, ParameterSemantic, auc_inf_is_reportable
from pharmacoweave.core.units import AMOUNT, DRUG_MASS, convert, parse_unit
from pharmacoweave.errors import UnsupportedDesign
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_profiles import ConcentrationProfile, TimedConcentration

# Keep this literal in sync with tests/numerical_benchmarks/reference_linlog.py.
_ADJUSTED_R2_TIE = 1e-6
_TIME_KEY_DIGITS = 9


@dataclass(frozen=True)
class _Regression:
    n_points: int
    lambda_z: float
    adjusted_r2: float


def calculate_profile(profile: ConcentrationProfile, config: NcaConfig) -> NcaParameterSet:
    """Calculate one NCA parameter set. Partial results stay partial."""

    if profile.route is not Route.IV_BOLUS:
        raise UnsupportedDesign(
            f"{profile.subject_id} route {profile.route.value} is outside the locked iv_bolus calculator."
        )
    warnings: list[str] = []
    if profile.unit_error:
        return _empty(
            profile,
            config,
            excluded_reason=profile.unit_error,
            warnings=warnings,
        )
    _note_blq(profile.points, warnings)
    post = [point for point in profile.points if point.time_h >= 0]
    duplicate = _duplicate_time(post)
    if duplicate:
        return _empty(
            profile,
            config,
            excluded_reason=(
                "Duplicate sample times were found. This calculator does not average technical replicates."
            ),
            warnings=warnings,
        )
    if any(point.status is ObservationStatus.ABOVE_ULOQ for point in post):
        return _empty(
            profile,
            config,
            excluded_reason="An above-ULOQ sample is present. This calculator does not impute the ULOQ.",
            warnings=warnings,
        )
    quantified = [point for point in post if point.status is ObservationStatus.QUANTIFIED and point.value is not None]
    quantified.sort(key=lambda point: point.time_h)
    predose_quantified = [
        point
        for point in profile.points
        if point.time_h < 0 and point.status is ObservationStatus.QUANTIFIED and point.value is not None
    ]
    if not quantified:
        reason = "No quantified concentrations."
        if predose_quantified:
            reason = "A quantified predose concentration is outside the single-dose IV bolus calculator."
        return _empty(profile, config, excluded_reason=reason, warnings=warnings)

    cmax_value = max(point.value or 0.0 for point in quantified)
    tmax_point = next(point for point in quantified if point.value == cmax_value)
    unit = profile.concentration_unit or tmax_point.value_unit or ""
    common = {
        "cmax": cmax_value,
        "cmax_unit": unit,
        "tmax": tmax_point.time_h,
        "tmax_unit": "h",
    }
    if predose_quantified:
        warnings.append(
            "A quantified predose concentration is outside the single-dose IV bolus calculator. "
            "AUC and the terminal phase were omitted."
        )
        return _record(profile, config, warnings=warnings, excluded_reason=warnings[-1], **common)

    internal = _internal_censor(post, quantified)
    if internal and config.blq_rule == "keep_as_censored":
        warnings.append(
            "A non-quantified sample lies between quantified samples. "
            "It was left censored and was not imputed, so AUC and the terminal phase were omitted."
        )
        return _record(profile, config, warnings=warnings, excluded_reason=warnings[-1], **common)
    if internal and config.blq_rule == "exclude":
        warnings.append(
            "A non-quantified sample between quantified samples was excluded. "
            "The trapezoid bridges that gap and can bias AUC."
        )

    auc_last, _c0 = _auc_last(quantified, warnings)
    regression = _terminal_regression(quantified, tmax_point.time_h, config.lambda_z_min_points)
    qualified = (
        regression is not None
        and regression.lambda_z > 0
        and regression.adjusted_r2 >= config.lambda_z_min_adjusted_r2
    )
    _note_terminal(regression, quantified, tmax_point.time_h, config, qualified, warnings)
    auc_last_unit = _auc_unit(unit) if auc_last is not None and unit else None
    terminal_n = regression.n_points if regression is not None else _points_after(quantified, tmax_point.time_h)
    adjusted = _storable_r2(regression)

    auc_inf = None
    auc_inf_unit = None
    fraction = None
    clearance = None
    clearance_unit = None
    clearance_semantic = None
    volume = None
    volume_unit = None
    volume_semantic = None
    if qualified and regression is not None and auc_last is not None and auc_last_unit is not None:
        last = quantified[-1].value or 0.0
        extrapolated = last / regression.lambda_z
        total = auc_last + extrapolated
        fraction_value = extrapolated / total if total > 0 else 1.0
        if fraction_value <= config.auc_inf_max_extrapolated_fraction:
            auc_inf = total
            auc_inf_unit = auc_last_unit
            fraction = fraction_value
            warnings.append(
                "AUC∞ extrapolation fraction was "
                f"{fraction_value:.6g}, within the reporting limit "
                f"{config.auc_inf_max_extrapolated_fraction:.6g}."
            )
            derived = _clearance_and_volume(
                dose=profile.dose_value,
                dose_unit=profile.dose_unit,
                auc=total,
                auc_unit=auc_last_unit,
                lambda_z=regression.lambda_z,
            )
            if derived is None:
                warnings.append(
                    "Clearance and Vz were omitted because the dose unit and the concentration unit "
                    "do not share a mass or amount dimension."
                )
            else:
                clearance, clearance_unit, volume, volume_unit = derived
                clearance_semantic = ParameterSemantic.CL
                volume_semantic = ParameterSemantic.VZ
        else:
            warnings.append(
                "AUC∞ extrapolation fraction was "
                f"{fraction_value:.6g}, above the reporting limit "
                f"{config.auc_inf_max_extrapolated_fraction:.6g}. "
                "AUC∞, clearance, and Vz were omitted. λz was kept because the terminal regression met its rule."
            )

    record = _record(
        profile,
        config,
        warnings=warnings,
        terminal_phase_qualified=qualified,
        auc_last=auc_last,
        auc_last_unit=auc_last_unit,
        auc_inf=auc_inf,
        auc_inf_unit=auc_inf_unit,
        auc_inf_extrap_fraction=fraction,
        lambda_z=regression.lambda_z if qualified and regression is not None else None,
        lambda_z_unit="1/h" if qualified and regression is not None else None,
        half_life=(math.log(2.0) / regression.lambda_z) if qualified and regression is not None else None,
        half_life_unit="h" if qualified and regression is not None else None,
        terminal_phase_n_points=terminal_n,
        terminal_phase_adjusted_r2=adjusted,
        clearance=clearance,
        clearance_unit=clearance_unit,
        clearance_semantic=clearance_semantic,
        volume=volume,
        volume_unit=volume_unit,
        volume_semantic=volume_semantic,
        **common,
    )
    if record.auc_inf is not None and not auc_inf_is_reportable(
        record, max_extrapolated_fraction=config.auc_inf_max_extrapolated_fraction
    ):
        raise RuntimeError("AUC∞ was stored even though the reporting rule rejected it.")
    return record


def _auc_last(quantified: list[TimedConcentration], warnings: list[str]) -> tuple[float | None, float | None]:
    if len(quantified) < 2:
        warnings.append("Fewer than two quantified points. AUC was omitted.")
        return None, None
    total = 0.0
    first, second = quantified[0], quantified[1]
    c0 = None
    if first.time_h > 0:
        backed = _back_extrapolate(first, second)
        if backed is None:
            warnings.append(
                "C0 was not back-extrapolated because the first interval does not decline. "
                "AUC starts at the first quantified sample."
            )
        else:
            c0, piece = backed
            total += piece
            warnings.append(
                f"Back-extrapolated C0 was {c0:.6g} {first.value_unit} for the interval from 0 to the first sample. "
                "C0 is not Cmax."
            )
    for left, right in zip(quantified, quantified[1:], strict=False):
        total += _trapezoid(left.value or 0.0, right.value or 0.0, right.time_h - left.time_h)
    return total, c0


def _trapezoid(c1: float, c2: float, dt: float) -> float:
    if dt <= 0:
        raise ValueError("NCA intervals require increasing time.")
    if c1 > 0 and c2 > 0 and c2 < c1:
        return (c1 - c2) * dt / math.log(c1 / c2)
    return 0.5 * (c1 + c2) * dt


def _back_extrapolate(first: TimedConcentration, second: TimedConcentration) -> tuple[float, float] | None:
    c1 = first.value or 0.0
    c2 = second.value or 0.0
    if not (c1 > 0 and c2 > 0 and c2 < c1 and second.time_h > first.time_h):
        return None
    slope = math.log(c2 / c1) / (second.time_h - first.time_h)
    if slope >= 0:
        return None
    c0 = c1 * math.exp(-slope * first.time_h)
    return c0, (c0 - c1) / (-slope)


def _terminal_regression(
    quantified: list[TimedConcentration],
    tmax_h: float,
    min_points: int,
) -> _Regression | None:
    after = [point for point in quantified if point.time_h > tmax_h and (point.value or 0.0) > 0]
    chosen: _Regression | None = None
    if len(after) < min_points:
        return None
    for size in range(min_points, len(after) + 1):
        window = after[-size:]
        fitted = _fit_window(window)
        if fitted is None:
            continue
        if chosen is None or _prefer(fitted, chosen):
            chosen = fitted
    return chosen


def _fit_window(window: list[TimedConcentration]) -> _Regression | None:
    n = len(window)
    if n < 3:
        return None
    times = [point.time_h for point in window]
    logs = [math.log(point.value or 0.0) for point in window]
    mean_t = sum(times) / n
    mean_y = sum(logs) / n
    ss_tt = sum((time - mean_t) ** 2 for time in times)
    ss_tot = sum((value - mean_y) ** 2 for value in logs)
    if ss_tt <= 0 or ss_tot <= 0:
        return None
    ss_ty = sum((time - mean_t) * (value - mean_y) for time, value in zip(times, logs, strict=True))
    slope = ss_ty / ss_tt
    if slope >= 0:
        return None
    intercept = mean_y - slope * mean_t
    ss_res = sum((value - (intercept + slope * time)) ** 2 for time, value in zip(times, logs, strict=True))
    r2 = 1.0 - ss_res / ss_tot
    if r2 > 1.0:
        r2 = 1.0
    adjusted = 1.0 - (1.0 - r2) * (n - 1) / (n - 2)
    return _Regression(n_points=n, lambda_z=-slope, adjusted_r2=adjusted)


def _prefer(candidate: _Regression, incumbent: _Regression) -> bool:
    if abs(candidate.adjusted_r2 - incumbent.adjusted_r2) <= _ADJUSTED_R2_TIE:
        return candidate.n_points > incumbent.n_points
    return candidate.adjusted_r2 > incumbent.adjusted_r2


def _note_terminal(
    regression: _Regression | None,
    quantified: list[TimedConcentration],
    tmax_h: float,
    config: NcaConfig,
    qualified: bool,
    warnings: list[str],
) -> None:
    available = _points_after(quantified, tmax_h)
    if regression is None:
        warnings.append(
            f"Terminal phase is not qualified. {available} quantified point(s) after Tmax "
            f"and no declining regression met the minimum of {config.lambda_z_min_points} points. "
            "λz, t½, and AUC∞ were omitted."
        )
        return
    state = "qualified" if qualified else "not qualified"
    warnings.append(
        f"Terminal phase is {state}. Regression used {regression.n_points} points after Tmax "
        f"({available} available). Adjusted R² was {regression.adjusted_r2:.6g}; "
        f"the minimum is {config.lambda_z_min_adjusted_r2:.6g}."
    )
    if not qualified:
        warnings.append("λz, t½, and AUC∞ were omitted because the terminal phase is not qualified.")


def _storable_r2(regression: _Regression | None) -> float | None:
    if regression is None:
        return None
    value = min(1.0, regression.adjusted_r2)
    if value < 0:
        return None
    return value


def _points_after(quantified: list[TimedConcentration], tmax_h: float) -> int:
    return sum(1 for point in quantified if point.time_h > tmax_h and (point.value or 0.0) > 0)


def _clearance_and_volume(
    *,
    dose: float,
    dose_unit: str,
    auc: float,
    auc_unit: str,
    lambda_z: float,
) -> tuple[float, str, float, str] | None:
    dose_dimensions = parse_unit(dose_unit).dimensions
    auc_dimensions = parse_unit(auc_unit).dimensions
    if dose_dimensions == DRUG_MASS and auc_dimensions == (1, 1, 0, -1, 0):
        dose_mg = convert(dose, dose_unit, "mg")
        auc_canon = convert(auc, auc_unit, "mg*h/L")
    elif dose_dimensions == AMOUNT and auc_dimensions == (1, 0, 0, -1, 1):
        dose_mg = convert(dose, dose_unit, "mol")
        auc_canon = convert(auc, auc_unit, "mol*h/L")
    else:
        return None
    if auc_canon <= 0 or lambda_z <= 0:
        return None
    clearance = dose_mg / auc_canon
    volume = clearance / lambda_z
    return clearance, "L/h", volume, "L"


def _note_blq(points: tuple[TimedConcentration, ...], warnings: list[str]) -> None:
    predose = sum(point.status is ObservationStatus.BLQ and point.time_h < 0 for point in points)
    post = sum(point.status is ObservationStatus.BLQ and point.time_h >= 0 for point in points)
    if predose:
        warnings.append(f"{predose} predose BLQ sample(s) were ignored.")
    if post:
        noun = "sample was" if post == 1 else "samples were"
        warnings.append(
            f"{post} post-dose BLQ {noun} left censored. "
            "They were not used as concentrations and were not set to zero."
        )


def _duplicate_time(points: list[TimedConcentration]) -> bool:
    seen: set[float] = set()
    for point in points:
        key = round(point.time_h, _TIME_KEY_DIGITS)
        if key in seen:
            return True
        seen.add(key)
    return False


def _internal_censor(post: list[TimedConcentration], quantified: list[TimedConcentration]) -> bool:
    if len(quantified) < 2:
        return False
    start = quantified[0].time_h
    end = quantified[-1].time_h
    return any(start < point.time_h < end and point.status is not ObservationStatus.QUANTIFIED for point in post)


def _auc_unit(concentration_unit: str) -> str:
    if "/" not in concentration_unit:
        return f"{concentration_unit}*h"
    left, right = concentration_unit.split("/", 1)
    return f"{left}*h/{right}"


def _empty(
    profile: ConcentrationProfile,
    config: NcaConfig,
    *,
    excluded_reason: str,
    warnings: list[str],
) -> NcaParameterSet:
    return _record(
        profile,
        config,
        warnings=warnings,
        excluded_reason=excluded_reason,
        terminal_phase_qualified=False,
    )


def _record(profile: ConcentrationProfile, config: NcaConfig, **fields: object) -> NcaParameterSet:
    payload = {
        "subject_id": profile.subject_id,
        "analyte_id": profile.analyte_id,
        "route": profile.route,
        "dose_event_id": profile.dose_event_id,
        "auc_method": config.auc_method,
        "terminal_phase_qualified": False,
    }
    payload.update(fields)
    return NcaParameterSet.model_validate(payload)
