"""Log-linear fit of a one-compartment IV bolus.

The regression is ln(C) = ln(Dose / V) - (CL / V) * t on quantified points.
It recovers CL and V when the curve is monoexponential. It does not fit
two-compartment models, absorption, or a separate additive error model.
BLQ points are omitted rather than replaced with zero.
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
from pharmacoweave.core.profiles import ConcentrationProfile, iv_bolus_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.specs import ONE_COMPARTMENT_IV_BOLUS
from pharmacoweave.version import __version__

_LIMITATIONS = (
    "Log-linear ordinary least squares assumes a one-compartment IV bolus and proportional residuals.",
    "BLQ concentrations were omitted and were not replaced with zero.",
    "Research use only. The estimates are not a dosing recommendation.",
)


@dataclass(frozen=True)
class IvBolusFit:
    """One subject-level fit and its diagnostic record."""

    fit: ModelFit
    estimates: tuple[ParameterEstimate, ...]
    diagnostic: DiagnosticRecord
    observation_ids: tuple[str, ...]
    predicted_mg_l: tuple[float, ...]


def fit_iv_bolus(study: PharmacologyStudy) -> list[IvBolusFit]:
    """Fit each single-dose IV bolus profile in ``study``."""

    return [fit_iv_bolus_profile(profile) for profile in iv_bolus_profiles(study)]


def fit_iv_bolus_profile(profile: ConcentrationProfile) -> IvBolusFit:
    """Fit one profile. A non-declining curve is a failed diagnostic, not a clearance."""

    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("fit_iv_bolus_profile requires route iv_bolus.")
    fit_id = f"fit-{profile.subject_id}-{profile.analyte_id}"
    dataset_id = f"{profile.subject_id}:{profile.analyte_id}:{profile.dose_event_id}"
    points = profile.points
    if len(points) < 3 or profile.dose_mg <= 0:
        diagnostic = DiagnosticRecord(
            fit_id=fit_id,
            status=DiagnosticStatus.FAIL,
            limitations=[*_LIMITATIONS, "Fewer than three quantified points, or the dose is not positive."],
            warnings=["CL and V were not estimated."],
        )
        fit = ModelFit(
            fit_id=fit_id,
            model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
            dataset_id=dataset_id,
            software_version=__version__,
            estimation_method="log_linear_ols",
            context_of_use=ONE_COMPARTMENT_IV_BOLUS.context_of_use,
            diagnostic_status=DiagnosticStatus.FAIL,
            converged=False,
            limitations=list(diagnostic.limitations),
        )
        return IvBolusFit(fit, (), diagnostic, (), ())

    xs = [point.time_h for point in points]
    ys = [math.log(point.concentration_mg_l) for point in points]
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    intercept = ybar - slope * xbar
    kel = -slope
    identifiable = kel > 0 and math.exp(intercept) > 0
    if not identifiable:
        diagnostic = DiagnosticRecord(
            fit_id=fit_id,
            status=DiagnosticStatus.FAIL,
            limitations=[*_LIMITATIONS, "The log-linear slope is not a declining concentration."],
            warnings=["CL and V were not estimated."],
        )
        fit = ModelFit(
            fit_id=fit_id,
            model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
            dataset_id=dataset_id,
            software_version=__version__,
            estimation_method="log_linear_ols",
            context_of_use=ONE_COMPARTMENT_IV_BOLUS.context_of_use,
            diagnostic_status=DiagnosticStatus.FAIL,
            converged=False,
            limitations=list(diagnostic.limitations),
        )
        return IvBolusFit(fit, (), diagnostic, (), ())

    volume = profile.dose_mg / math.exp(intercept)
    clearance = kel * volume
    predicted = tuple((profile.dose_mg / volume) * math.exp(-(clearance / volume) * point.time_h) for point in points)
    relative = [
        abs(pred - point.concentration_mg_l) / point.concentration_mg_l
        for pred, point in zip(predicted, points, strict=True)
    ]
    max_relative = max(relative)
    sst = sum((y - ybar) ** 2 for y in ys)
    sse = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys, strict=True))
    r_squared = 1.0 if sst <= 1e-18 else 1.0 - max(0.0, sse) / sst
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(_LIMITATIONS),
        warnings=[] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."],
        metrics=[
            Metric(name="log_r_squared", value=r_squared),
            Metric(name="max_relative_residual", value=max_relative),
            Metric(name="n_points", value=float(n)),
        ],
    )
    fit = ModelFit(
        fit_id=fit_id,
        model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=dataset_id,
        software_version=__version__,
        estimation_method="log_linear_ols",
        context_of_use=ONE_COMPARTMENT_IV_BOLUS.context_of_use,
        diagnostic_status=status,
        converged=True,
        limitations=list(_LIMITATIONS),
    )
    estimates = (
        ParameterEstimate(
            fit_id=fit_id,
            parameter="CL",
            semantic=ParameterSemantic.CL,
            value=clearance,
            unit="L/h",
            estimation_method="log_linear_ols",
            route=Route.IV_BOLUS,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V",
            semantic=ParameterSemantic.V,
            value=volume,
            unit="L",
            estimation_method="log_linear_ols",
            route=Route.IV_BOLUS,
        ),
    )
    return IvBolusFit(
        fit=fit,
        estimates=estimates,
        diagnostic=diagnostic,
        observation_ids=tuple(point.observation_id for point in points),
        predicted_mg_l=predicted,
    )
