"""Weighted nonlinear least squares for a one-compartment constant-rate IV infusion.

The mean model is the closed form in ``iv_infusion_concentration``. CL and V are
estimated by Levenberg–Marquardt in log(CL) and log(V). Weights are
``uniform``, ``1/y``, or ``1/y^2``. The solver does not integrate an ODE.

A missing infusion duration, a duration that is not shorter than the last
quantified time, or a curve that does not identify both parameters is
``DiagnosticStatus.FAIL`` and has no clearance. BLQ points are omitted rather
than replaced with zero. Estimates are CL and V, not CL/F or V/F.
"""

import math
from dataclasses import dataclass

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import (
    DiagnosticRecord,
    DiagnosticStatus,
    Metric,
    ModelFit,
    ParameterEstimate,
    ParameterSemantic,
)
from pharmacoweave.core.numerics import student_t_two_sided_quantile
from pharmacoweave.core.profiles import ConcentrationProfile, iv_infusion_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import iv_infusion_concentration
from pharmacoweave.pk_models.fit_iv import ParameterUncertainty
from pharmacoweave.pk_models.specs import ONE_COMPARTMENT_IV_INFUSION
from pharmacoweave.version import __version__

_WEIGHTS = {
    "uniform": "weighted_nls:uniform",
    "1/y": "weighted_nls:1/y",
    "1/y^2": "weighted_nls:1/y^2",
}
_CI_LEVEL = 0.95
_MIN_CURVATURE_RATIO = 1e-8
_LIMITATIONS = (
    "Weighted nonlinear least squares assumes a single-dose one-compartment constant-rate IV infusion.",
    "Elimination is linear. The mean model is iv_infusion_concentration, not a numerical ODE.",
    "CL and V are absolute intravenous parameters, not CL/F or V/F.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton log-parameter approximation.",
    "Research use only. The estimates are not a dosing recommendation.",
)


@dataclass(frozen=True)
class InfusionFit:
    """One subject-level constant-rate IV infusion fit and its diagnostic record.

    ``predicted_mg_l`` is empty when the fit fails. Failed fits do not carry
    a clearance. Times and observations are the quantified points only.
    """

    fit: ModelFit
    estimates: tuple[ParameterEstimate, ...]
    diagnostic: DiagnosticRecord
    observation_ids: tuple[str, ...]
    predicted_mg_l: tuple[float, ...]
    times_h: tuple[float, ...]
    observed_mg_l: tuple[float, ...]
    subject_id: str


def fit_one_compartment_iv_infusion(
    study: PharmacologyStudy,
    *,
    weight: str = "1/y^2",
) -> list[InfusionFit]:
    """Fit each single-dose constant-rate IV infusion profile."""

    return [
        fit_one_compartment_iv_infusion_profile(profile, weight=weight) for profile in iv_infusion_profiles(study)
    ]


def fit_one_compartment_iv_infusion_profile(
    profile: ConcentrationProfile,
    *,
    weight: str = "1/y^2",
) -> InfusionFit:
    """Fit one infusion profile. Duration and identifiability failures omit CL and V."""

    method = _method(weight)
    if profile.route is not Route.IV_INFUSION:
        raise AnalysisError("fit_one_compartment_iv_infusion_profile requires route iv_infusion.")
    limitations = (*_LIMITATIONS, f"Weights are {weight} on quantified concentrations.")
    if profile.interval_regimen != "single_dose" or profile.n_doses != 1:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="This fit accepts a single-dose constant-rate IV infusion profile.",
        )
    duration = profile.infusion_duration_h
    if duration is None or duration <= 0 or not math.isfinite(duration):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Infusion duration is missing or not positive.",
        )
    if len(profile.points) < 3 or profile.dose_mg <= 0:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than three quantified points, or the dose is not positive.",
        )
    times = [point.time_h for point in profile.points]
    observations = [point.concentration_mg_l for point in profile.points]
    if any(observation <= 0 or not math.isfinite(observation) for observation in observations):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="A quantified concentration is not positive.",
        )
    if max(times) - min(times) <= 1e-12:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Quantified sample times do not vary, so CL and V are not identifiable.",
        )
    if not (max(times) > duration):
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Infusion duration is not shorter than the observation window.",
        )
    if max(observations) / min(observations) < 1.01:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="Quantified concentrations do not show enough curvature to identify CL and V.",
        )
    slope_reason, kel0, intercept0 = _post_infusion_line(times, observations, duration)
    if slope_reason is not None:
        return _failed(profile, method=method, limitations=limitations, reason=slope_reason)

    fitted = _fit_parameters(
        times,
        observations,
        dose=profile.dose_mg,
        duration=duration,
        weight=weight,
        kel=kel0,
        intercept=intercept0,
    )
    if fitted is None:
        return _failed(
            profile,
            method=method,
            limitations=limitations,
            reason="The weighted nonlinear least squares fit did not converge, or CL and V are not identifiable.",
        )
    clearance, volume, weighted_sse, iterations, inv_cl, inv_v = fitted
    predicted = tuple(
        _reference_concentration(
            dose=profile.dose_mg,
            duration=duration,
            clearance=clearance,
            volume=volume,
            time_h=time_h,
        )
        for time_h in times
    )
    relative = [
        abs(pred - observed) / observed for pred, observed in zip(predicted, observations, strict=True) if observed > 0
    ]
    max_relative = max(relative) if relative else math.inf
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
        limitations=(*limitations, f"The infusion duration used in the mean model is {duration:.6g} h."),
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
            Metric(name="infusion_duration_h", value=duration, unit="h"),
        ],
        clearance_uncertainty=cl_uncertainty,
        volume_uncertainty=v_uncertainty,
    )


def _reference_concentration(*, dose: float, duration: float, clearance: float, volume: float, time_h: float) -> float:
    """Evaluate the public infusion closed form in mg/L."""

    return iv_infusion_concentration(
        dose=dose,
        dose_unit="mg",
        infusion_duration=duration,
        infusion_duration_unit="h",
        volume=volume,
        volume_unit="L",
        clearance=clearance,
        clearance_unit="L/h",
        time=time_h,
        time_unit="h",
        concentration_unit="mg/L",
    )


def _concentration_and_log_columns(
    dose: float,
    duration: float,
    clearance: float,
    volume: float,
    time_h: float,
) -> tuple[float, float, float] | None:
    """Return C and derivatives with respect to log(CL) and log(V).

    The concentration matches ``iv_infusion_concentration``. Columns are
    CL * dC/dCL and V * dC/dV for the log-parameter step.
    """

    if dose < 0 or duration <= 0 or clearance <= 0 or volume <= 0 or time_h < 0:
        return None
    rate = dose / duration
    kel = clearance / volume
    if dose == 0 or time_h == 0:
        return 0.0, 0.0, 0.0
    if time_h <= duration:
        decay = math.exp(-kel * time_h)
        concentration = -(rate / clearance) * math.expm1(-kel * time_h)
        d_clearance = -concentration / clearance + (rate / clearance) * decay * (time_h / volume)
        d_volume = -(rate * time_h / (volume * volume)) * decay
    else:
        end_decay = math.exp(-kel * duration)
        concentration_end = -(rate / clearance) * math.expm1(-kel * duration)
        tau = time_h - duration
        post = math.exp(-kel * tau)
        concentration = concentration_end * post
        d_end_clearance = -concentration_end / clearance + (rate / clearance) * end_decay * (duration / volume)
        d_end_volume = -(rate * duration / (volume * volume)) * end_decay
        d_clearance = d_end_clearance * post - concentration * tau / volume
        d_volume = d_end_volume * post + concentration * tau * clearance / (volume * volume)
    column_cl = d_clearance * clearance
    column_v = d_volume * volume
    if not (math.isfinite(concentration) and math.isfinite(column_cl) and math.isfinite(column_v)):
        return None
    return concentration, column_cl, column_v


def _post_infusion_line(
    times: list[float],
    observations: list[float],
    duration: float,
) -> tuple[str | None, float | None, float | None]:
    """Screen points after the infusion stops.

    Two or more post-infusion points must decline. Fewer than two points do
    not reject the curve; the optimizer then has to identify CL and V from
    the remaining samples.
    """

    post = [(time_h, observation) for time_h, observation in zip(times, observations, strict=True) if time_h > duration]
    if len(post) < 2:
        return None, None, None
    xs = [time_h for time_h, _observation in post]
    ys = [math.log(observation) for _time_h, observation in post]
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0 or not math.isfinite(sxx):
        return "Post-infusion sample times do not vary, so CL and V are not identifiable.", None, None
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    intercept = ybar - slope * xbar
    kel = -slope
    if not (math.isfinite(kel) and math.isfinite(intercept) and kel > 0):
        return "Post-infusion concentrations are not declining, so CL and V are not identifiable.", None, None
    return None, kel, intercept


def _fit_parameters(
    times: list[float],
    observations: list[float],
    *,
    dose: float,
    duration: float,
    weight: str,
    kel: float | None,
    intercept: float | None,
) -> tuple[float, float, float, int, float | None, float | None] | None:
    starts: list[tuple[float, float]] = []
    if kel is not None and intercept is not None:
        analytical = _analytical_start(dose, duration, kel, intercept)
        if analytical is not None:
            starts.append(analytical)
    ranked = sorted(
        _grid_starts(times, observations, dose=dose, duration=duration, weight=weight),
        key=_sse_of_start(times, observations, dose=dose, duration=duration, weight=weight),
    )
    for candidate in ranked[:3]:
        if all(not _same_start(candidate, existing) for existing in starts):
            starts.append(candidate)
    best: tuple[float, float, float, int, float | None, float | None] | None = None
    for clearance0, volume0 in starts:
        fitted = _levenberg_marquardt(
            times,
            observations,
            dose=dose,
            duration=duration,
            clearance=clearance0,
            volume=volume0,
            weight=weight,
        )
        if fitted is None:
            continue
        if best is None or fitted[2] < best[2]:
            best = fitted
    return best


def _analytical_start(dose: float, duration: float, kel: float, intercept: float) -> tuple[float, float] | None:
    c_end = math.exp(intercept - kel * duration)
    factor = -math.expm1(-kel * duration)
    if not (c_end > 0 and factor > 0 and math.isfinite(c_end) and math.isfinite(factor)):
        return None
    clearance = (dose / duration) * factor / c_end
    volume = clearance / kel
    if not (clearance > 0 and volume > 0 and math.isfinite(clearance) and math.isfinite(volume)):
        return None
    return clearance, volume


def _grid_starts(
    times: list[float],
    observations: list[float],
    *,
    dose: float,
    duration: float,
    weight: str,
) -> list[tuple[float, float]]:
    rate = dose / duration
    cl_scale = rate / max(observations)
    starts: list[tuple[float, float]] = []
    for i_cl in range(11):
        clearance = cl_scale * 10.0 ** (-2.0 + 4.0 * i_cl / 10.0)
        for i_kel in range(11):
            kel = 10.0 ** (-3.0 + 4.0 * i_kel / 10.0)
            volume = clearance / kel
            trial = _sse(times, observations, dose, duration, clearance, volume, weight)
            if _finite_positive(clearance, volume) and trial < math.inf:
                starts.append((clearance, volume))
    return starts


def _sse_of_start(
    times: list[float],
    observations: list[float],
    *,
    dose: float,
    duration: float,
    weight: str,
):
    def key(start: tuple[float, float]) -> float:
        return _sse(times, observations, dose, duration, start[0], start[1], weight)

    return key


def _same_start(left: tuple[float, float], right: tuple[float, float]) -> bool:
    return abs(left[0] - right[0]) <= 1e-8 * max(1.0, abs(right[0])) and abs(left[1] - right[1]) <= 1e-8 * max(
        1.0, abs(right[1])
    )


def _finite_positive(*values: float) -> bool:
    return all(math.isfinite(value) and value > 0 for value in values)


def _weight_of(weight: str, observation: float) -> float:
    if weight == "uniform":
        return 1.0
    if weight == "1/y":
        return 1.0 / observation
    return 1.0 / (observation * observation)


def _sse(
    times: list[float],
    observations: list[float],
    dose: float,
    duration: float,
    clearance: float,
    volume: float,
    weight: str,
) -> float:
    total = 0.0
    for time_h, observation in zip(times, observations, strict=True):
        fitted = _concentration_and_log_columns(dose, duration, clearance, volume, time_h)
        if fitted is None:
            return math.inf
        residual = observation - fitted[0]
        total += _weight_of(weight, observation) * residual * residual
    return total if math.isfinite(total) else math.inf


def _levenberg_marquardt(
    times: list[float],
    observations: list[float],
    *,
    dose: float,
    duration: float,
    clearance: float,
    volume: float,
    weight: str,
) -> tuple[float, float, float, int, float | None, float | None] | None:
    """Return CL, V, weighted SSE, iterations, and the inverse diagonal of JᵀWJ."""

    def sse(clearance_value: float, volume_value: float) -> float:
        return _sse(times, observations, dose, duration, clearance_value, volume_value, weight)

    def information(clearance_value: float, volume_value: float) -> tuple[float, float, float, float, float] | None:
        a = b = c = grad_cl = grad_v = 0.0
        for time_h, observation in zip(times, observations, strict=True):
            fitted = _concentration_and_log_columns(dose, duration, clearance_value, volume_value, time_h)
            if fitted is None:
                return None
            concentration, column_cl, column_v = fitted
            residual = observation - concentration
            point_weight = _weight_of(weight, observation)
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
        info = information(clearance, volume)
        if info is None:
            return None
        a, b, c, grad_cl, grad_v = info
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
            if not (math.isfinite(step_cl) and math.isfinite(step_v)):
                lam = min(lam * 8.0, 1e12)
                continue
            if max(abs(step_cl), abs(step_v)) <= 1e-10:
                converged = True
                accepted = True
                break
            if max(abs(step_cl), abs(step_v)) > 5.0:
                lam = min(lam * 8.0, 1e12)
                continue
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
    info_final = information(clearance, volume)
    if info_final is None:
        return None
    a_final, b_final, c_final, _g_cl, _g_v = info_final
    if not converged or _curvature_ratio(a_final, b_final, c_final) < _MIN_CURVATURE_RATIO:
        return None
    positive = clearance > 0 and volume > 0
    finite = math.isfinite(clearance) and math.isfinite(volume) and math.isfinite(current)
    if not (positive and finite):
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
    if larger <= 1e-30 or smaller < 0.0:
        return 0.0
    return smaller / larger


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
    return (
        _one_log_uncertainty(clearance, math.sqrt(sigma2 * inv_cl), t_crit),
        _one_log_uncertainty(volume, math.sqrt(sigma2 * inv_v), t_crit),
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


def _succeeded(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    status: DiagnosticStatus,
    clearance: float,
    volume: float,
    predicted: tuple[float, ...],
    warnings: list[str],
    metrics: list[Metric],
    clearance_uncertainty: ParameterUncertainty,
    volume_uncertainty: ParameterUncertainty,
) -> InfusionFit:
    fit_id = _fit_id(profile)
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(limitations),
        warnings=warnings,
        metrics=metrics,
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_IV_INFUSION.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_IV_INFUSION.context_of_use,
        diagnostic_status=status,
        converged=True,
        limitations=list(limitations),
    )
    estimates = (
        ParameterEstimate(
            fit_id=fit_id,
            parameter="CL",
            semantic=ParameterSemantic.CL,
            value=clearance,
            unit="L/h",
            estimation_method=method,
            route=Route.IV_INFUSION,
            standard_error=clearance_uncertainty.standard_error,
            ci_lower=clearance_uncertainty.ci_lower,
            ci_upper=clearance_uncertainty.ci_upper,
            ci_level=clearance_uncertainty.ci_level,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V",
            semantic=ParameterSemantic.V,
            value=volume,
            unit="L",
            estimation_method=method,
            route=Route.IV_INFUSION,
            standard_error=volume_uncertainty.standard_error,
            ci_lower=volume_uncertainty.ci_lower,
            ci_upper=volume_uncertainty.ci_upper,
            ci_level=volume_uncertainty.ci_level,
        ),
    )
    return _pack(profile, fit, estimates, diagnostic, predicted)


def _failed(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    reason: str,
) -> InfusionFit:
    fit_id = _fit_id(profile)
    stated = [*limitations, reason]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=["CL and V were not estimated."],
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_IV_INFUSION.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_IV_INFUSION.context_of_use,
        diagnostic_status=DiagnosticStatus.FAIL,
        converged=False,
        limitations=stated,
    )
    return _pack(profile, fit, (), diagnostic, ())


def _pack(
    profile: ConcentrationProfile,
    fit: ModelFit,
    estimates: tuple[ParameterEstimate, ...],
    diagnostic: DiagnosticRecord,
    predicted: tuple[float, ...],
) -> InfusionFit:
    return InfusionFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        observation_ids=tuple(point.observation_id for point in profile.points),
        predicted_mg_l=predicted,
        times_h=tuple(point.time_h for point in profile.points),
        observed_mg_l=tuple(point.concentration_mg_l for point in profile.points),
        subject_id=profile.subject_id,
    )


def _fit_id(profile: ConcentrationProfile) -> str:
    return f"fit-{profile.subject_id}-{profile.analyte_id}"


def _dataset_id(profile: ConcentrationProfile) -> str:
    return f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}"
