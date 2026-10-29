"""Independent linear-up/log-down reference.

This transcription is the external numeric reference for the locked IV bolus
NCA. It does not import PharmacoWeave, PKNCA, or Phoenix. The trapezoid is the
linear-up/log-down rule those tools document. The terminal slope is the
maximum adjusted R² on points strictly after Tmax. Adjusted R² values within
1e-6 are tied, and the tie goes to more points. Keep that 1e-6 literal aligned
with ``pharmacoweave.pk_nca.calculator``.
"""

import math

_ADJUSTED_R2_TIE = 1e-6


def reference_nca(times_h: list[float], concentrations: list[float], *, min_points: int = 3) -> dict[str, float | int]:
    """Return trapezoidal NCA for quantified points only. BLQ is not an input."""

    if len(times_h) != len(concentrations) or len(times_h) < 2:
        raise ValueError("The reference needs at least two quantified points.")
    order = sorted(range(len(times_h)), key=times_h.__getitem__)
    times = [times_h[index] for index in order]
    concs = [concentrations[index] for index in order]
    auc_last = 0.0
    c0 = concs[0]
    if times[0] > 0 and concs[0] > concs[1] > 0:
        slope = math.log(concs[1] / concs[0]) / (times[1] - times[0])
        c0 = concs[0] * math.exp(-slope * times[0])
        auc_last += (c0 - concs[0]) / (-slope)
    for index in range(len(times) - 1):
        auc_last += _trapezoid(concs[index], concs[index + 1], times[index + 1] - times[index])
    cmax = max(concs)
    tmax = times[concs.index(cmax)]
    after_times = [time for time, conc in zip(times, concs, strict=True) if time > tmax]
    after_concs = [conc for time, conc in zip(times, concs, strict=True) if time > tmax]
    chosen = _best_terminal(after_times, after_concs, min_points)
    if chosen is None:
        raise ValueError("The reference terminal phase was not qualified.")
    n_points, lambda_z, adjusted_r2 = chosen
    extrapolated = concs[-1] / lambda_z
    auc_inf = auc_last + extrapolated
    return {
        "c0": c0,
        "auc_last": auc_last,
        "auc_inf": auc_inf,
        "fraction": extrapolated / auc_inf,
        "cmax": cmax,
        "tmax_h": tmax,
        "lambda_z": lambda_z,
        "half_life_h": math.log(2.0) / lambda_z,
        "n_points": n_points,
        "adjusted_r2": adjusted_r2,
    }


def _trapezoid(c1: float, c2: float, dt: float) -> float:
    if c1 > 0 and c2 > 0 and c2 < c1:
        return (c1 - c2) * dt / math.log(c1 / c2)
    return 0.5 * (c1 + c2) * dt


def _best_terminal(
    times: list[float],
    concs: list[float],
    min_points: int,
) -> tuple[int, float, float] | None:
    chosen: tuple[int, float, float] | None = None
    if len(times) < min_points:
        return None
    for size in range(min_points, len(times) + 1):
        fitted = _fit(times[-size:], concs[-size:])
        if fitted is None:
            continue
        if chosen is None or _prefer(fitted, chosen):
            chosen = fitted
    return chosen


def _fit(times: list[float], concs: list[float]) -> tuple[int, float, float] | None:
    n = len(times)
    logs = [math.log(value) for value in concs]
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
    r2 = min(1.0, 1.0 - ss_res / ss_tot)
    adjusted = 1.0 - (1.0 - r2) * (n - 1) / (n - 2)
    return n, -slope, adjusted


def _prefer(candidate: tuple[int, float, float], incumbent: tuple[int, float, float]) -> bool:
    if abs(candidate[2] - incumbent[2]) <= _ADJUSTED_R2_TIE:
        return candidate[0] > incumbent[0]
    return candidate[2] > incumbent[2]
