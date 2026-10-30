"""Weighted nonlinear least squares for a one-compartment IV bolus.

The mean model is C(t) = (Dose / V) * exp(-(CL / V) * t). Weights are
``1/y`` or ``1/y^2`` on quantified concentrations. The solver is
Levenberg–Marquardt in log(CL) and log(V), so a successful step stays
positive. It does not call an external optimizer.

A curve that is not declining on the log-linear screen, or that does not
identify both parameters, is ``DiagnosticStatus.FAIL`` and has no clearance.
BLQ points are omitted rather than replaced with zero.
"""

import math

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import DiagnosticStatus, Metric
from pharmacoweave.core.numerics import student_t_two_sided_quantile
from pharmacoweave.core.profiles import ConcentrationProfile, iv_bolus_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.fit_iv import IvBolusFit, ParameterUncertainty, _failed, _succeeded, screen_log_linear

_WEIGHTS = {"1/y": "weighted_nls:1/y", "1/y^2": "weighted_nls:1/y^2"}
_CI_LEVEL = 0.95
_LIMITATIONS = (
    "Weighted nonlinear least squares assumes a one-compartment IV bolus.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton log-parameter approximation.",
    "Research use only. The estimates are not a dosing recommendation.",
)


def fit_iv_bolus_wnls(study: PharmacologyStudy, *, weight: str = "1/y^2") -> list[IvBolusFit]:
    """Fit each single-dose IV bolus profile with weighted nonlinear least squares."""

    return [fit_iv_bolus_profile_wnls(profile, weight=weight) for profile in iv_bolus_profiles(study)]


def fit_iv_bolus_profile_wnls(profile: ConcentrationProfile, *, weight: str = "1/y^2") -> IvBolusFit:
    """Fit one profile. Non-declining and non-identifiable curves do not return CL or V."""

    method = _method(weight)
    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("fit_iv_bolus_profile_wnls requires route iv_bolus.")
    limitations = (*_LIMITATIONS, f"Weights are {weight} on quantified concentrations.")
    if len(profile.points) < 3 or profile.dose_mg <= 0:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than three quantified points, or the dose is not positive.",
        )
    screen = screen_log_linear(profile.points)
    if not screen.identifiable or not screen.declining or screen.kel is None or screen.intercept is None:
        return _failed(profile, method=method, limitations=limitations, reason=screen.reason)

    volume0 = profile.dose_mg / math.exp(screen.intercept)
    clearance0 = screen.kel * volume0
    times = [point.time_h for point in profile.points]
    observations = [point.concentration_mg_l for point in profile.points]
    fitted = _levenberg_marquardt(
        times,
        observations,
        dose=profile.dose_mg,
        clearance=clearance0,
        volume=volume0,
        weight=weight,
    )
    if fitted is None:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="The weighted nonlinear least squares fit did not converge, or CL and V are not identifiable.",
        )
    clearance, volume, weighted_sse, iterations, inv_cl, inv_v = fitted
    predicted = tuple(_concentration(profile.dose_mg, clearance, volume, time_h) for time_h in times)
    relative = [
        abs(pred - observed) / observed for pred, observed in zip(predicted, observations, strict=True)
    ]
    max_relative = max(relative)
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."]
    df = len(times) - 2
    cl_uncertainty, v_uncertainty = _log_parameter_uncertainties(
        clearance=clearance,
        volume=volume,
        weighted_sse=weighted_sse,
        df=df,
        inv_cl=inv_cl,
        inv_v=inv_v,
        warnings=warnings,
    )
    return _succeeded(
        profile,
        method=method,
        limitations=limitations,
        status=status,
        clearance=clearance,
        volume=volume,
        predicted=predicted,
        warnings=warnings,
        metrics=[
            Metric(name="weighted_sse", value=weighted_sse),
            Metric(name="max_relative_residual", value=max_relative),
            Metric(name="n_points", value=float(len(times))),
            Metric(name="n_iterations", value=float(iterations)),
        ],
        clearance_uncertainty=cl_uncertainty,
        volume_uncertainty=v_uncertainty,
    )


def _log_parameter_uncertainties(
    *,
    clearance: float,
    volume: float,
    weighted_sse: float,
    df: int,
    inv_cl: float | None,
    inv_v: float | None,
    warnings: list[str],
) -> tuple[ParameterUncertainty, ParameterUncertainty]:
    if df < 1 or inv_cl is None or inv_v is None or inv_cl < 0 or inv_v < 0:
        warnings.append("Information matrix or residual degrees of freedom did not support confidence intervals.")
        return ParameterUncertainty(), ParameterUncertainty()
    sigma2 = max(0.0, weighted_sse) / float(df)
    if not math.isfinite(sigma2):
        warnings.append("Residual variance was not finite; confidence intervals were omitted.")
        return ParameterUncertainty(), ParameterUncertainty()
    t_crit = student_t_two_sided_quantile(df, _CI_LEVEL)
    se_log_cl = math.sqrt(sigma2 * inv_cl)
    se_log_v = math.sqrt(sigma2 * inv_v)
    return (
        _one_log_uncertainty(clearance, se_log_cl, t_crit),
        _one_log_uncertainty(volume, se_log_v, t_crit),
    )


def _one_log_uncertainty(value: float, se_log: float, t_crit: float) -> ParameterUncertainty:
    half_width = min(t_crit * se_log, 50.0)
    return ParameterUncertainty(
        standard_error=value * se_log,
        ci_lower=value * math.exp(-half_width),
        ci_upper=value * math.exp(half_width),
        ci_level=_CI_LEVEL,
    )


def _method(weight: str) -> str:
    try:
        return _WEIGHTS[weight]
    except KeyError as exc:
        allowed = ", ".join(sorted(_WEIGHTS))
        raise AnalysisError(f"weight must be one of {allowed}.") from exc


def _concentration(dose: float, clearance: float, volume: float, time_h: float) -> float:
    return (dose / volume) * math.exp(-(clearance / volume) * time_h)


def _levenberg_marquardt(
    times: list[float],
    observations: list[float],
    *,
    dose: float,
    clearance: float,
    volume: float,
    weight: str,
) -> tuple[float, float, float, int, float | None, float | None] | None:
    """Return CL, V, weighted SSE, iterations, and inverse diagonal of JᵀWJ."""

    def weight_of(observation: float) -> float:
        if weight == "1/y":
            return 1.0 / observation
        return 1.0 / (observation * observation)

    def sse(clearance_value: float, volume_value: float) -> float:
        total = 0.0
        for time_h, observation in zip(times, observations, strict=True):
            residual = observation - _concentration(dose, clearance_value, volume_value, time_h)
            total += weight_of(observation) * residual * residual
        return total

    def information(clearance_value: float, volume_value: float) -> tuple[float, float, float, float, float]:
        a = b = c = grad_cl = grad_v = 0.0
        for time_h, observation in zip(times, observations, strict=True):
            fitted = _concentration(dose, clearance_value, volume_value, time_h)
            d_clearance = -fitted * time_h / volume_value
            d_volume = fitted * (clearance_value * time_h - volume_value) / (volume_value * volume_value)
            column_cl = d_clearance * clearance_value
            column_v = d_volume * volume_value
            residual = observation - fitted
            point_weight = weight_of(observation)
            a += point_weight * column_cl * column_cl
            b += point_weight * column_cl * column_v
            c += point_weight * column_v * column_v
            grad_cl += point_weight * column_cl * residual
            grad_v += point_weight * column_v * residual
        return a, b, c, grad_cl, grad_v

    lam = 1e-3
    current = sse(clearance, volume)
    if not math.isfinite(current):
        return None
    converged = False
    iterations = 0
    for _step in range(50):
        iterations += 1
        a, b, c, grad_cl, grad_v = information(clearance, volume)
        if max(abs(grad_cl), abs(grad_v)) <= 1e-8 * (1.0 + current):
            converged = True
            break
        accepted = False
        for _attempt in range(20):
            scale_a = a if a > 0.0 else 1.0
            scale_c = c if c > 0.0 else 1.0
            damped_a = a + lam * scale_a
            damped_c = c + lam * scale_c
            determinant = damped_a * damped_c - b * b
            if abs(determinant) <= 1e-30 * (1.0 + abs(damped_a) * abs(damped_c)):
                lam = min(lam * 8.0, 1e12)
                continue
            step_cl = (grad_cl * damped_c - b * grad_v) / determinant
            step_v = (damped_a * grad_v - b * grad_cl) / determinant
            if max(abs(step_cl), abs(step_v)) <= 1e-10:
                converged = True
                accepted = True
                break
            updated_cl = clearance * math.exp(step_cl)
            updated_v = volume * math.exp(step_v)
            if not (math.isfinite(updated_cl) and math.isfinite(updated_v)):
                lam = min(lam * 8.0, 1e12)
                continue
            trial = sse(updated_cl, updated_v)
            if not math.isfinite(trial):
                lam = min(lam * 8.0, 1e12)
                continue
            if trial <= current * (1.0 + 1e-14):
                clearance, volume, current = updated_cl, updated_v, trial
                lam = max(lam / 3.0, 1e-14)
                accepted = True
                break
            lam = min(lam * 8.0, 1e12)
        if converged or not accepted:
            break
    a_final, b_final, c_final, _g_cl, _g_v = information(clearance, volume)
    curvature = _curvature_ratio(a_final, b_final, c_final)
    if not converged or curvature < 1e-8 or not (clearance > 0 and volume > 0):
        return None
    if not (math.isfinite(clearance) and math.isfinite(volume) and math.isfinite(current)):
        return None
    det = a_final * c_final - b_final * b_final
    inv_cl: float | None = None
    inv_v: float | None = None
    if det > 1e-30 * (1.0 + abs(a_final * c_final)):
        inv_cl = c_final / det
        inv_v = a_final / det
    return clearance, volume, current, iterations, inv_cl, inv_v


def _curvature_ratio(a: float, b: float, c: float) -> float:
    """Return the smaller eigenvalue of JᵀWJ divided by the larger one."""

    trace = a + c
    determinant = a * c - b * b
    discriminant = trace * trace - 4.0 * determinant
    if discriminant < 0.0:
        discriminant = 0.0
    root = math.sqrt(discriminant)
    larger = 0.5 * (trace + root)
    smaller = 0.5 * (trace - root)
    if larger <= 1e-30:
        return 0.0
    return smaller / larger
