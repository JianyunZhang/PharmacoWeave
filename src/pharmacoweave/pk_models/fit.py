"""One-compartment IV bolus estimation of CL and V.

The structural model is ``C(t) = (Dose / V) * exp(-(CL / V) * t)``. Parameters
are fit by nonlinear least squares on the logarithms of CL and V. A fit that
does not converge, hits a bound, or is not identifiable is recorded as a
failure and does not carry parameter estimates.

The objective uses quantified post-dose concentrations only. BLQ samples are
omitted; this is not a censored likelihood. Research use only.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from pharmacoweave.core.analysis_export import PredictionRow
from pharmacoweave.core.entities import ObservationStatus, Route
from pharmacoweave.core.modeling import (
    AnalysisRecords,
    DiagnosticRecord,
    DiagnosticStatus,
    ErrorModel,
    Metric,
    ModelFit,
    ParameterEstimate,
    ParameterSemantic,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import AMOUNT, CONC_AMOUNT, CONC_MASS, DRUG_MASS, convert, parse_unit
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import ValidationFailed
from pharmacoweave.pk_models.specs import ONE_COMPARTMENT_IV_BOLUS
from pharmacoweave.pk_profiles import ConcentrationProfile, TimedConcentration, load_iv_bolus_profiles
from pharmacoweave.version import __version__

_CL_MIN = 1e-8
_CL_MAX = 1e6
_V_MIN = 1e-8
_V_MAX = 1e6
_WALD_Z = 1.959963984540054
_CONTEXT = "Individual one-compartment IV bolus fit for research use. Not a dosing recommendation."
_LIMITS = (
    "Research use only. This fit is not a dosing recommendation, a medical-device output, or a regulatory validation.",
    "The structural model is one-compartment IV bolus with linear elimination. "
    "Goodness of fit does not prove that mechanism.",
    "BLQ and other non-quantified samples were omitted from the objective. This is not a censored-data likelihood.",
)


@dataclass(frozen=True)
class FittedCurve:
    """Predicted concentrations on a grid, in the profile concentration unit."""

    subject_id: str
    analyte_id: str
    dose_event_id: str
    points: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class IvBolusFit:
    """A fit attempt. Failed profiles contribute diagnostics and no estimates."""

    records: AnalysisRecords
    predictions: tuple[PredictionRow, ...]
    curves: tuple[FittedCurve, ...]


def fit_one_compartment_iv_bolus(
    study: PharmacologyStudy,
    *,
    error_model: str = "proportional",
    max_nfev: int = 200,
    initial: tuple[float, float] | None = None,
) -> IvBolusFit:
    """Estimate absolute CL (L/h) and V (L) for each IV bolus profile.

    ``initial`` is an optional ``(CL, V)`` start in L/h and L. It is not a
    prior. ``max_nfev`` limits the optimizer. Non-convergence and
    non-identifiability omit estimates and cannot be marked pass.
    """

    if error_model not in {ErrorModel.PROPORTIONAL.value, ErrorModel.ADDITIVE.value}:
        raise ValueError("error_model must be 'proportional' or 'additive'.")
    if max_nfev < 1:
        raise ValueError("max_nfev must be at least 1.")
    if initial is not None and (initial[0] <= 0 or initial[1] <= 0):
        raise ValueError("initial CL and V must be positive.")
    report = validate_study(study)
    if not report.ok:
        raise ValidationFailed(report)
    profiles = load_iv_bolus_profiles(study)
    method = f"nonlinear least squares, {error_model} error, log-parameter parameterization"
    fits: list[ModelFit] = []
    diagnostics: list[DiagnosticRecord] = []
    estimates: list[ParameterEstimate] = []
    predictions: list[PredictionRow] = []
    curves: list[FittedCurve] = []
    for profile in profiles:
        outcome = _fit_profile(
            profile,
            study_id=study.metadata.study_id,
            error_model=error_model,
            method=method,
            max_nfev=max_nfev,
            initial=initial,
        )
        fits.append(outcome.fit)
        diagnostics.append(outcome.diagnostic)
        estimates.extend(outcome.estimates)
        predictions.extend(outcome.predictions)
        if outcome.curve is not None:
            curves.append(outcome.curve)
    records = AnalysisRecords(
        model_specs=[ONE_COMPARTMENT_IV_BOLUS],
        model_fits=fits,
        parameter_estimates=estimates,
        diagnostics=diagnostics,
    )
    return IvBolusFit(records=records, predictions=tuple(predictions), curves=tuple(curves))


@dataclass(frozen=True)
class _Outcome:
    fit: ModelFit
    diagnostic: DiagnosticRecord
    estimates: tuple[ParameterEstimate, ...]
    predictions: tuple[PredictionRow, ...]
    curve: FittedCurve | None


def _fit_profile(
    profile: ConcentrationProfile,
    *,
    study_id: str,
    error_model: str,
    method: str,
    max_nfev: int,
    initial: tuple[float, float] | None,
) -> _Outcome:
    fit_id = f"FIT-{profile.subject_id}-{profile.analyte_id}"
    warnings = _sample_warnings(profile.points)
    scale = _scale(profile)
    quantified = [
        point
        for point in profile.points
        if point.time_h >= 0 and point.status is ObservationStatus.QUANTIFIED and point.value is not None
    ]
    failure = _precheck(profile, quantified, scale)
    if failure is not None:
        return _fail(fit_id, study_id, method, failure, warnings, converged=None)

    assert scale is not None
    dose, canon_unit = scale
    assert profile.concentration_unit is not None
    times = np.array([point.time_h for point in quantified], dtype=float)
    observed = np.array(
        [convert(point.value or 0.0, profile.concentration_unit, canon_unit) for point in quantified],
        dtype=float,
    )
    if initial is None:
        start = _default_start(times, observed, dose)
    else:
        start = initial
    start = (_interior(start[0], _CL_MIN, _CL_MAX), _interior(start[1], _V_MIN, _V_MAX))
    proportional = error_model == ErrorModel.PROPORTIONAL.value

    def residual(x: np.ndarray) -> np.ndarray:
        cl = math.exp(float(x[0]))
        volume = math.exp(float(x[1]))
        predicted = _predict(dose, cl, volume, times)
        if proportional:
            return (observed - predicted) / predicted
        return observed - predicted

    try:
        result = least_squares(
            residual,
            np.array([math.log(start[0]), math.log(start[1])], dtype=float),
            bounds=([math.log(_CL_MIN), math.log(_V_MIN)], [math.log(_CL_MAX), math.log(_V_MAX)]),
            max_nfev=max_nfev,
        )
    except (ValueError, FloatingPointError, RuntimeError) as exc:
        return _fail(
            fit_id,
            study_id,
            method,
            f"The optimizer raised {type(exc).__name__}. Estimates were omitted.",
            [*warnings, str(exc)],
            converged=False,
        )

    clearance = math.exp(float(result.x[0]))
    volume = math.exp(float(result.x[1]))
    jacobian = np.asarray(result.jac, dtype=float)
    identifiable = _identifiable(jacobian)
    on_bound = bool(np.any(np.asarray(result.active_mask) != 0))
    hit_bound = on_bound or _near_bound(clearance, _CL_MIN, _CL_MAX) or _near_bound(volume, _V_MIN, _V_MAX)
    if not bool(result.success) or not math.isfinite(clearance) or not math.isfinite(volume):
        return _fail(
            fit_id,
            study_id,
            method,
            "The optimizer did not converge. Estimates were omitted.",
            warnings,
            converged=False,
        )
    if hit_bound:
        return _fail(
            fit_id,
            study_id,
            method,
            "A parameter hit its bound. Estimates were omitted and the fit was not marked pass.",
            warnings,
            converged=True,
        )
    if not identifiable:
        return _fail(
            fit_id,
            study_id,
            method,
            "CL and V are not identifiable from this profile. Estimates were omitted.",
            warnings,
            converged=True,
        )

    objective = np.asarray(result.fun, dtype=float)
    rss = float(np.sum(np.square(objective)))
    standard_errors = _standard_errors(jacobian, rss, clearance, volume)
    limitations = list(_LIMITS)
    if standard_errors is None:
        limitations.append(
            "Standard errors were not estimated because the residual degrees of freedom "
            "were zero or the local curvature was not invertible."
        )
    else:
        limitations.append(
            "Standard errors and Wald intervals are a local quadratic approximation on the "
            "log-parameter scale, converted with the delta method. They are not a regulatory interval."
        )
    estimates = _estimates(
        fit_id,
        method,
        clearance,
        volume,
        standard_errors,
        warnings,
    )
    prediction_rows, curve = _predictions(
        profile,
        quantified,
        fit_id=fit_id,
        dose=dose,
        clearance=clearance,
        volume=volume,
        canon_unit=canon_unit,
        proportional=proportional,
    )
    metrics = [
        Metric(name="n_quantified", value=float(len(quantified))),
        Metric(name="weighted_sum_of_squares", value=rss),
    ]
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=study_id,
        software_version=__version__,
        estimation_method=method,
        context_of_use=_CONTEXT,
        diagnostic_status=DiagnosticStatus.PASS,
        converged=True,
        limitations=limitations,
    )
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.PASS,
        limitations=limitations,
        warnings=warnings,
        metrics=metrics,
    )
    return _Outcome(fit, diagnostic, estimates, prediction_rows, curve)


def _precheck(
    profile: ConcentrationProfile,
    quantified: list[TimedConcentration],
    scale: tuple[float, str] | None,
) -> str | None:
    if profile.unit_error:
        return f"{profile.unit_error} Estimates were omitted."
    if profile.route is not Route.IV_BOLUS:
        return "The one-compartment IV bolus fit does not estimate CL for a non-bolus route."
    if scale is None:
        return "Dose and concentration dimensions do not match, so CL and V were not estimated."
    if len(quantified) < 2:
        return "Fewer than two quantified post-dose concentrations. CL and V are not identifiable."
    span = max(point.time_h for point in quantified) - min(point.time_h for point in quantified)
    if span <= 1e-9:
        return "Quantified samples do not span time. CL and V are not identifiable."
    keys = {round(point.time_h, 9) for point in quantified}
    if len(keys) < len(quantified):
        return "Duplicate sample times were found. Replicates are not averaged, and the fit was omitted."
    return None


def _scale(profile: ConcentrationProfile) -> tuple[float, str] | None:
    if profile.concentration_unit is None or profile.unit_error is not None:
        return None
    try:
        dose_dimensions = parse_unit(profile.dose_unit).dimensions
        conc_dimensions = parse_unit(profile.concentration_unit).dimensions
    except ValueError:
        return None
    if dose_dimensions == DRUG_MASS and conc_dimensions == CONC_MASS:
        return convert(profile.dose_value, profile.dose_unit, "mg"), "mg/L"
    if dose_dimensions == AMOUNT and conc_dimensions == CONC_AMOUNT:
        return convert(profile.dose_value, profile.dose_unit, "mol"), "mol/L"
    return None


def _default_start(times: np.ndarray, concs: np.ndarray, dose: float) -> tuple[float, float]:
    if float(times[-1]) > float(times[0]) and float(concs[-1]) < float(concs[0]) and float(concs[0]) > 0:
        kel = -(math.log(float(concs[-1])) - math.log(float(concs[0]))) / (float(times[-1]) - float(times[0]))
        c0 = float(concs[0]) * math.exp(kel * float(times[0]))
        if kel > 0 and c0 > 0:
            return kel * dose / c0, dose / c0
    volume = dose / float(concs[0])
    return 0.1 * volume, volume


def _interior(value: float, lower: float, upper: float) -> float:
    return min(max(value, lower * 1.0000001), upper / 1.0000001)


def _near_bound(value: float, lower: float, upper: float, *, relative: float = 1e-3) -> bool:
    """Treat a parameter pressed against an optimizer bound as a bound hit.

    ``least_squares`` can stop a hair inside the bound when the objective is
    already flat, leaving ``active_mask`` clear. A clearance or volume that
    close to the artificial bound is not an interior estimate.
    """

    return value <= lower * (1.0 + relative) or value >= upper / (1.0 + relative)


def _predict(dose: float, clearance: float, volume: float, times: np.ndarray) -> np.ndarray:
    predicted = (dose / volume) * np.exp(-(clearance / volume) * times)
    return np.maximum(predicted, 1e-30)


def _identifiable(jacobian: np.ndarray) -> bool:
    if jacobian.ndim != 2 or jacobian.shape[1] != 2 or jacobian.shape[0] < 2:
        return False
    singular = np.linalg.svd(jacobian, compute_uv=False)
    if singular.size < 2 or float(singular[0]) <= 0:
        return False
    return float(singular[-1]) / float(singular[0]) >= 1e-8


def _standard_errors(
    jacobian: np.ndarray,
    rss: float,
    clearance: float,
    volume: float,
) -> tuple[float, float] | None:
    dof = int(jacobian.shape[0]) - 2
    if dof <= 0:
        return None
    try:
        covariance = (rss / dof) * np.linalg.inv(jacobian.T @ jacobian)
    except np.linalg.LinAlgError:
        return None
    variance = np.diag(covariance)
    if not np.all(np.isfinite(variance)) or np.any(variance < 0):
        return None
    se_log = np.sqrt(variance)
    return float(clearance * se_log[0]), float(volume * se_log[1])


def _estimates(
    fit_id: str,
    method: str,
    clearance: float,
    volume: float,
    standard_errors: tuple[float, float] | None,
    warnings: list[str],
) -> tuple[ParameterEstimate, ...]:
    rows = []
    for name, semantic, value, se in (
        ("CL", ParameterSemantic.CL, clearance, None if standard_errors is None else standard_errors[0]),
        ("V", ParameterSemantic.V, volume, None if standard_errors is None else standard_errors[1]),
    ):
        lower = upper = level = None
        if se is not None and math.isfinite(se) and se >= 0:
            lower_value = value - _WALD_Z * se
            upper_value = value + _WALD_Z * se
            if lower_value > 0:
                lower = lower_value
                upper = upper_value
                level = 0.95
            else:
                warnings.append(f"The Wald interval for {name} crossed zero and was omitted.")
        rows.append(
            ParameterEstimate(
                fit_id=fit_id,
                parameter=name,
                semantic=semantic,
                value=value,
                unit="L/h" if semantic is ParameterSemantic.CL else "L",
                estimation_method=method,
                route=Route.IV_BOLUS,
                standard_error=se,
                ci_lower=lower,
                ci_upper=upper,
                ci_level=level,
                bound_hit=False,
            )
        )
    return tuple(rows)


def _predictions(
    profile: ConcentrationProfile,
    quantified: list[TimedConcentration],
    *,
    fit_id: str,
    dose: float,
    clearance: float,
    volume: float,
    canon_unit: str,
    proportional: bool,
) -> tuple[tuple[PredictionRow, ...], FittedCurve]:
    assert profile.concentration_unit is not None
    rows: list[PredictionRow] = []
    for point in quantified:
        observed_canon = convert(point.value or 0.0, profile.concentration_unit, canon_unit)
        predicted_canon = float(_predict(dose, clearance, volume, np.array([point.time_h], dtype=float))[0])
        predicted = convert(predicted_canon, canon_unit, profile.concentration_unit)
        observed = point.value or 0.0
        residual = observed - predicted
        if proportional:
            weighted = (observed_canon - predicted_canon) / predicted_canon
        else:
            weighted = residual
        rows.append(
            PredictionRow(
                fit_id=fit_id,
                observation_id=point.observation_id,
                subject_id=profile.subject_id,
                time_h=point.time_h,
                observed=observed,
                predicted=predicted,
                residual=residual,
                weighted_residual=weighted,
                unit=profile.concentration_unit,
            )
        )
    t_last = max(point.time_h for point in quantified)
    grid = tuple(t_last * index / 79 for index in range(80))
    curve_values = _predict(dose, clearance, volume, np.array(grid, dtype=float))
    curve_points = tuple(
        (time_h, convert(float(value), canon_unit, profile.concentration_unit))
        for time_h, value in zip(grid, curve_values, strict=True)
    )
    curve = FittedCurve(
        subject_id=profile.subject_id,
        analyte_id=profile.analyte_id,
        dose_event_id=profile.dose_event_id,
        points=curve_points,
    )
    return tuple(rows), curve


def _sample_warnings(points: tuple[TimedConcentration, ...]) -> list[str]:
    post_blq = sum(point.status is ObservationStatus.BLQ and point.time_h >= 0 for point in points)
    predose = sum(point.time_h < 0 and point.status is ObservationStatus.QUANTIFIED for point in points)
    warnings: list[str] = []
    if post_blq:
        warnings.append(
            f"{post_blq} post-dose BLQ sample(s) were omitted from the objective and were not set to zero."
        )
    if predose:
        warnings.append("Quantified predose concentrations were omitted from the objective.")
    return warnings


def _fail(
    fit_id: str,
    study_id: str,
    method: str,
    reason: str,
    warnings: list[str],
    *,
    converged: bool | None,
) -> _Outcome:
    limitations = [reason, *_LIMITS]
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=study_id,
        software_version=__version__,
        estimation_method=method,
        context_of_use=_CONTEXT,
        diagnostic_status=DiagnosticStatus.FAIL,
        converged=converged,
        limitations=limitations,
    )
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=limitations,
        warnings=warnings,
    )
    return _Outcome(fit, diagnostic, (), (), None)
