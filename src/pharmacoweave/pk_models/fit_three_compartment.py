"""Noise-free recovery of a three-compartment IV bolus.

Parameters are CL, V1, Q2, V2, Q3, and V3. A failed recovery leaves every
clearance and volume empty.
"""

import math
from dataclasses import dataclass
from typing import Literal

from pharmacoweave.core.entities import Record, Route
from pharmacoweave.core.modeling import (
    DiagnosticRecord,
    DiagnosticStatus,
    Metric,
    ModelFit,
    ParameterEstimate,
    ParameterSemantic,
)
from pharmacoweave.core.numerics import (
    invert_symmetric,
    normalized_curvature_ratio,
    solve_linear_system,
    student_t_two_sided_quantile,
)
from pharmacoweave.core.profiles import ConcentrationProfile, iv_bolus_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import three_compartment_bolus_mg_l
from pharmacoweave.pk_models.fit_iv import ParameterUncertainty
from pharmacoweave.pk_models.specs import THREE_COMPARTMENT_IV_BOLUS
from pharmacoweave.version import __version__

_MIN_POINTS = 8
_LOG_STEP = 1e-6
_CI_LEVEL = 0.95
_MIN_CURVATURE = 1e-8
_WEIGHTS = {"uniform": "weighted_nls:uniform", "1/y": "weighted_nls:1/y", "1/y^2": "weighted_nls:1/y^2"}
_FIT_LIMITATIONS = (
    "Assumes a single-dose three-compartment IV bolus with central elimination.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Standard errors and 95% confidence intervals use the local Gauss-Newton log-parameter approximation.",
    "Research use only. The estimates are not a dosing recommendation.",
)


@dataclass(frozen=True)
class ThreeCompartmentFit:
    """One subject-level three-compartment IV bolus fit."""

    fit: ModelFit
    estimates: tuple[ParameterEstimate, ...]
    diagnostic: DiagnosticRecord
    observation_ids: tuple[str, ...]
    predicted_mg_l: tuple[float, ...]
    times_h: tuple[float, ...]
    observed_mg_l: tuple[float, ...]
    subject_id: str


class ThreeCompartmentRecovery(Record):
    """Recovered three-compartment parameters. Empty values mean the fit failed."""

    dose_mg: float
    clearance_l_h: float | None = None
    central_volume_l: float | None = None
    q2_l_h: float | None = None
    v2_l: float | None = None
    q3_l_h: float | None = None
    v3_l: float | None = None
    converged: bool = False
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def recover_three_compartment_iv_bolus(
    *,
    dose_mg: float,
    times_h: list[float],
    concentrations_mg_l: list[float],
) -> ThreeCompartmentRecovery:
    """Recover parameters from a noise-free central concentration curve."""

    limitations = [
        "Recovery uses the three-compartment IV bolus closed form.",
        "A failed recovery does not report clearance or volume.",
        "The result is not a dosing recommendation.",
    ]
    if (
        not math.isfinite(dose_mg)
        or dose_mg <= 0
        or len(times_h) != len(concentrations_mg_l)
        or len(times_h) < _MIN_POINTS
        or any(not math.isfinite(time) or time < 0 for time in times_h)
        or any(later <= earlier for earlier, later in zip(times_h, times_h[1:], strict=False))
        or any(not math.isfinite(value) or value <= 0 for value in concentrations_mg_l)
        or times_h[0] != 0.0
    ):
        return _failed(
            dose_mg,
            limitations,
            "The curve needs a positive dose and at least eight increasing positive points, starting at time zero.",
        )
    if concentrations_mg_l[-1] >= concentrations_mg_l[0]:
        return _failed(dose_mg, limitations, "The curve does not decline, so no clearance is reported.")
    best: tuple[float, tuple[float, float, float, float, float, float]] | None = None
    for seed in _seeds(dose_mg, times_h, concentrations_mg_l):
        fitted = _fit_seed(seed, times_h, concentrations_mg_l, dose_mg, weights=None)
        if fitted is None:
            continue
        sse, params = fitted
        if best is None or sse < best[0]:
            best = (sse, params)
    scale = math.fsum(value * value for value in concentrations_mg_l)
    if best is None or best[0] > 1e-8 * max(scale, 1.0):
        return _failed(
            dose_mg,
            limitations,
            "The closed form did not recover a finite parameter set, so no clearance is reported.",
        )
    cl, v1, q2, v2, q3, v3 = best[1]
    return ThreeCompartmentRecovery(
        dose_mg=dose_mg,
        clearance_l_h=cl,
        central_volume_l=v1,
        q2_l_h=q2,
        v2_l=v2,
        q3_l_h=q3,
        v3_l=v3,
        converged=True,
        limitations=limitations,
    )


def _failed(dose_mg: float, limitations: list[str], reason: str) -> ThreeCompartmentRecovery:
    safe_dose = dose_mg if math.isfinite(dose_mg) and dose_mg > 0 else 0.0
    return ThreeCompartmentRecovery(dose_mg=safe_dose, limitations=[*limitations, reason])


def _seeds(
    dose_mg: float,
    times: list[float],
    concentrations: list[float],
) -> list[tuple[float, float, float, float, float, float]]:
    terminal = _terminal_slope(times, concentrations)
    auc = _auc_with_tail(times, concentrations, terminal)
    clearance = dose_mg / auc if auc > 0 else concentrations[0]
    volume = dose_mg / concentrations[0]
    seeds: list[tuple[float, float, float, float, float, float]] = []
    for q2_factor, v2_factor, q3_factor, v3_factor in (
        (1.0, 2.0, 0.3, 8.0),
        (0.5, 1.5, 0.2, 5.0),
        (2.0, 4.0, 0.5, 12.0),
        (1.5, 3.0, 0.1, 10.0),
        (0.8, 6.0, 0.4, 3.0),
    ):
        seeds.append(
            (
                clearance,
                volume,
                q2_factor * clearance,
                v2_factor * volume,
                q3_factor * clearance,
                v3_factor * volume,
            )
        )
    return seeds


def _terminal_slope(times: list[float], concentrations: list[float]) -> float:
    start = max(1, len(times) // 2)
    later_time, earlier_time = times[-1], times[start]
    later, earlier = concentrations[-1], concentrations[start]
    if later_time <= earlier_time or later >= earlier:
        return 0.1
    return max(1e-4, (math.log(earlier) - math.log(later)) / (later_time - earlier_time))


def _auc_with_tail(times: list[float], concentrations: list[float], slope: float) -> float:
    area = math.fsum(
        0.5 * (left + right) * (times[index + 1] - times[index])
        for index, (left, right) in enumerate(zip(concentrations, concentrations[1:], strict=False))
    )
    return area + concentrations[-1] / slope


def fit_three_compartment_iv_bolus(
    study: PharmacologyStudy,
    *,
    weight: str = "1/y^2",
) -> list[ThreeCompartmentFit]:
    """Fit each single-dose IV bolus profile with the three-compartment closed form."""

    return [fit_three_compartment_iv_bolus_profile(profile, weight=weight) for profile in iv_bolus_profiles(study)]


def fit_three_compartment_iv_bolus_profile(
    profile: ConcentrationProfile,
    *,
    weight: str = "1/y^2",
) -> ThreeCompartmentFit:
    """Fit one IV bolus profile. A failed fit returns no clearance."""

    method = _method(weight)
    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("fit_three_compartment_iv_bolus_profile requires route iv_bolus.")
    limitations = (*_FIT_LIMITATIONS, f"Weights are {weight} on quantified concentrations.")
    points = profile.points
    times = [point.time_h for point in points]
    observations = [point.concentration_mg_l for point in points]
    if len(points) < _MIN_POINTS or profile.dose_mg <= 0 or observations[-1] >= observations[0]:
        return _failed_fit(
            profile,
            method=method,
            limitations=limitations,
            reason="Fewer than eight quantified points, a nonpositive dose, or a curve that does not decline.",
        )
    weights = _observation_weights(weight, observations)
    best: tuple[float, tuple[float, float, float, float, float, float]] | None = None
    for seed in _seeds(profile.dose_mg, times, observations):
        fitted = _fit_seed(seed, times, observations, profile.dose_mg, weights=weights)
        if fitted is None:
            continue
        if best is None or fitted[0] < best[0]:
            best = fitted
    if best is None:
        return _failed_fit(
            profile,
            method=method,
            limitations=limitations,
            reason="The three-compartment fit did not converge, so no clearance is reported.",
        )
    sse, params = best
    logs = [math.log(value) for value in params]
    hessian = _normal_matrix(logs, times, observations, profile.dose_mg, weights)
    inverse = None if hessian is None else invert_symmetric(hessian)
    curvature = 0.0 if hessian is None else normalized_curvature_ratio(hessian)
    if inverse is None or curvature < _MIN_CURVATURE:
        return _failed_fit(
            profile,
            method=method,
            limitations=limitations,
            reason="The six parameters are not separately identifiable, so no clearance is reported.",
        )
    predicted = tuple(
        three_compartment_bolus_mg_l(
            dose_mg=profile.dose_mg,
            cl_l_h=params[0],
            v1_l=params[1],
            q2_l_h=params[2],
            v2_l=params[3],
            q3_l_h=params[4],
            v3_l=params[5],
            time_h=time,
        )
        for time in times
    )
    max_relative = max(abs(pred - obs) / obs for pred, obs in zip(predicted, observations, strict=True))
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."]
    uncertainties = _uncertainties(params, weighted_sse=sse, df=len(times) - 6, inverse=inverse, warnings=warnings)
    return _succeeded_fit(
        profile,
        method=method,
        limitations=limitations,
        status=status,
        params=params,
        uncertainties=uncertainties,
        predicted=predicted,
        warnings=warnings,
        metrics=[
            Metric(name="weighted_sse", value=sse),
            Metric(name="max_relative_residual", value=max_relative),
            Metric(name="curvature_ratio", value=curvature),
            Metric(name="n_points", value=float(len(times))),
        ],
    )


def _fit_seed(
    seed: tuple[float, float, float, float, float, float],
    times: list[float],
    observations: list[float],
    dose_mg: float,
    *,
    weights: list[float] | None,
) -> tuple[float, tuple[float, float, float, float, float, float]] | None:
    point_weights = weights if weights is not None else [1.0 for _ in observations]
    logs = [math.log(value) for value in seed]
    current = _sse(logs, times, observations, dose_mg, point_weights)
    if not math.isfinite(current):
        return None
    lam = 1e-3
    for _step in range(60):
        jacobian, residual = _jacobian(logs, times, observations, dose_mg, point_weights)
        if jacobian is None or residual is None:
            return None
        gradient = [
            math.fsum(row[index] * value for row, value in zip(jacobian, residual, strict=True))
            for index in range(6)
        ]
        normal = [
            [
                math.fsum(row[left] * row[right] for row in jacobian)
                for right in range(6)
            ]
            for left in range(6)
        ]
        if max(abs(value) for value in gradient) <= 1e-10 * (1.0 + current):
            break
        accepted = False
        for _attempt in range(16):
            damped = [
                [
                    value + (lam * max(normal[index][index], 1e-8) if left == index else 0.0)
                    for index, value in enumerate(row)
                ]
                for left, row in enumerate(normal)
            ]
            step = solve_linear_system(damped, gradient)
            if step is None or any(not math.isfinite(value) or abs(value) > 2.0 for value in step):
                lam = min(lam * 8.0, 1e12)
                continue
            trial_logs = [value + delta for value, delta in zip(logs, step, strict=True)]
            trial = _sse(trial_logs, times, observations, dose_mg, point_weights)
            if math.isfinite(trial) and trial <= current * (1.0 + 1e-12):
                logs, current = trial_logs, trial
                lam = max(lam / 3.0, 1e-14)
                accepted = True
                break
            lam = min(lam * 8.0, 1e12)
        if not accepted:
            break
    params = tuple(math.exp(value) for value in logs)
    if not all(math.isfinite(value) and value > 0 for value in params):
        return None
    return current, params  # type: ignore[return-value]


def _predict(logs: list[float], times: list[float], dose_mg: float) -> list[float] | None:
    cl, v1, q2, v2, q3, v3 = (math.exp(value) for value in logs)
    try:
        return [
            three_compartment_bolus_mg_l(
                dose_mg=dose_mg,
                cl_l_h=cl,
                v1_l=v1,
                q2_l_h=q2,
                v2_l=v2,
                q3_l_h=q3,
                v3_l=v3,
                time_h=time,
            )
            for time in times
        ]
    except ValueError:
        return None


def _sse(
    logs: list[float],
    times: list[float],
    observations: list[float],
    dose_mg: float,
    weights: list[float],
) -> float:
    predicted = _predict(logs, times, dose_mg)
    if predicted is None:
        return math.inf
    return math.fsum(
        weight * (obs - pred) ** 2
        for weight, obs, pred in zip(weights, observations, predicted, strict=True)
    )


def _jacobian(
    logs: list[float],
    times: list[float],
    observations: list[float],
    dose_mg: float,
    weights: list[float],
) -> tuple[list[list[float]] | None, list[float] | None]:
    base = _predict(logs, times, dose_mg)
    if base is None:
        return None, None
    columns: list[list[float]] = []
    for index in range(6):
        plus = logs.copy()
        minus = logs.copy()
        plus[index] += _LOG_STEP
        minus[index] -= _LOG_STEP
        high = _predict(plus, times, dose_mg)
        low = _predict(minus, times, dose_mg)
        if high is None or low is None:
            return None, None
        columns.append([(up - down) / (2.0 * _LOG_STEP) for up, down in zip(high, low, strict=True)])
    jacobian = [[columns[param][point] for param in range(6)] for point in range(len(times))]
    residual = [obs - pred for obs, pred in zip(observations, base, strict=True)]
    scales = [math.sqrt(weight) for weight in weights]
    jacobian = [[value * scales[point] for value in row] for point, row in enumerate(jacobian)]
    residual = [value * scales[point] for point, value in enumerate(residual)]
    return jacobian, residual


def _normal_matrix(
    logs: list[float],
    times: list[float],
    observations: list[float],
    dose_mg: float,
    weights: list[float],
) -> list[list[float]] | None:
    jacobian, _residual = _jacobian(logs, times, observations, dose_mg, weights)
    if jacobian is None:
        return None
    return [[math.fsum(row[left] * row[right] for row in jacobian) for right in range(6)] for left in range(6)]


def _observation_weights(weight: str, observations: list[float]) -> list[float]:
    if weight == "uniform":
        return [1.0 for _ in observations]
    if weight == "1/y":
        return [1.0 / observation for observation in observations]
    if weight == "1/y^2":
        return [1.0 / (observation * observation) for observation in observations]
    raise AnalysisError("weight must be uniform, 1/y, or 1/y^2.")


def _method(weight: str) -> str:
    try:
        return _WEIGHTS[weight]
    except KeyError as exc:
        raise AnalysisError("weight must be uniform, 1/y, or 1/y^2.") from exc


def _uncertainties(
    values: tuple[float, float, float, float, float, float],
    *,
    weighted_sse: float,
    df: int,
    inverse: list[list[float]],
    warnings: list[str],
) -> tuple[ParameterUncertainty, ...]:
    if df < 1:
        warnings.append("Residual degrees of freedom did not support confidence intervals.")
        return tuple(ParameterUncertainty() for _ in values)
    sigma2 = max(0.0, weighted_sse) / float(df)
    if not math.isfinite(sigma2):
        warnings.append("Residual variance was not finite; confidence intervals were omitted.")
        return tuple(ParameterUncertainty() for _ in values)
    critical = student_t_two_sided_quantile(df, _CI_LEVEL)
    found: list[ParameterUncertainty] = []
    for index, value in enumerate(values):
        se_log = math.sqrt(max(0.0, sigma2 * inverse[index][index]))
        half_width = min(critical * se_log, 50.0)
        found.append(
            ParameterUncertainty(
                standard_error=value * se_log,
                ci_lower=value * math.exp(-half_width),
                ci_upper=value * math.exp(half_width),
                ci_level=_CI_LEVEL,
            )
        )
    return tuple(found)


def _succeeded_fit(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    status: DiagnosticStatus,
    params: tuple[float, float, float, float, float, float],
    uncertainties: tuple[ParameterUncertainty, ...],
    predicted: tuple[float, ...],
    warnings: list[str],
    metrics: list[Metric],
) -> ThreeCompartmentFit:
    fit_id = f"fit-3c-iv-{profile.subject_id}-{profile.analyte_id}"
    names = ("CL", "V1", "Q2", "V2", "Q3", "V3")
    semantics = (
        ParameterSemantic.CL,
        ParameterSemantic.V1,
        ParameterSemantic.Q2,
        ParameterSemantic.V2,
        ParameterSemantic.Q3,
        ParameterSemantic.V3,
    )
    units = ("L/h", "L", "L/h", "L", "L/h", "L")
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(limitations),
        warnings=warnings,
        metrics=metrics,
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=THREE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}",
        software_version=__version__,
        estimation_method=method,
        context_of_use=THREE_COMPARTMENT_IV_BOLUS.context_of_use,
        diagnostic_status=status,
        converged=True,
        limitations=list(limitations),
    )
    estimates = tuple(
        ParameterEstimate(
            fit_id=fit_id,
            parameter=name,
            semantic=semantic,
            value=value,
            unit=unit,
            estimation_method=method,
            route=Route.IV_BOLUS,
            standard_error=uncertainty.standard_error,
            ci_lower=uncertainty.ci_lower,
            ci_upper=uncertainty.ci_upper,
            ci_level=uncertainty.ci_level,
        )
        for name, semantic, value, unit, uncertainty in zip(
            names, semantics, params, units, uncertainties, strict=True
        )
    )
    return _pack_fit(profile, fit, estimates, diagnostic, predicted)


def _failed_fit(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    reason: str,
) -> ThreeCompartmentFit:
    fit_id = f"fit-3c-iv-{profile.subject_id}-{profile.analyte_id}"
    stated = [*limitations, reason]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=["CL, V1, Q2, V2, Q3, and V3 were not estimated."],
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=THREE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}",
        software_version=__version__,
        estimation_method=method,
        context_of_use=THREE_COMPARTMENT_IV_BOLUS.context_of_use,
        diagnostic_status=DiagnosticStatus.FAIL,
        converged=False,
        limitations=stated,
    )
    return _pack_fit(profile, fit, (), diagnostic, ())


def _pack_fit(
    profile: ConcentrationProfile,
    fit: ModelFit,
    estimates: tuple[ParameterEstimate, ...],
    diagnostic: DiagnosticRecord,
    predicted: tuple[float, ...],
) -> ThreeCompartmentFit:
    return ThreeCompartmentFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        observation_ids=tuple(point.observation_id for point in profile.points),
        predicted_mg_l=predicted,
        times_h=tuple(point.time_h for point in profile.points),
        observed_mg_l=tuple(point.concentration_mg_l for point in profile.points),
        subject_id=profile.subject_id,
    )
