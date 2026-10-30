"""Log-linear fit of a one-compartment IV bolus.

The regression is ln(C) = ln(Dose / V) - (CL / V) * t on quantified points.
It recovers CL and V when the curve is monoexponential. It does not fit
two-compartment models, absorption, or a separate additive error model.
BLQ points are omitted rather than replaced with zero.

Weighted nonlinear least squares for the same equation lives in
``pharmacoweave.pk_models.wnls``.
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
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint, iv_bolus_profiles
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
class LogLinearScreen:
    """Whether ln(concentration) versus time can identify a declining bolus curve."""

    identifiable: bool
    declining: bool
    kel: float | None
    intercept: float | None
    slope: float | None
    r_squared: float | None
    reason: str


@dataclass(frozen=True)
class IvBolusFit:
    """One subject-level fit and its diagnostic record.

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


def fit_iv_bolus(study: PharmacologyStudy) -> list[IvBolusFit]:
    """Fit each single-dose IV bolus profile in ``study`` by log-linear OLS."""

    return [fit_iv_bolus_profile(profile) for profile in iv_bolus_profiles(study)]


def fit_iv_bolus_profile(profile: ConcentrationProfile) -> IvBolusFit:
    """Fit one profile. A non-declining curve is a failed diagnostic, not a clearance."""

    if profile.route is not Route.IV_BOLUS:
        raise AnalysisError("fit_iv_bolus_profile requires route iv_bolus.")
    if len(profile.points) < 3 or profile.dose_mg <= 0:
        return _failed(
            profile,
            method="log_linear_ols",
            limitations=_LIMITATIONS,
            reason="Fewer than three quantified points, or the dose is not positive.",
        )
    screen = screen_log_linear(profile.points)
    if not screen.identifiable or not screen.declining or screen.kel is None or screen.intercept is None:
        return _failed(
            profile,
            method="log_linear_ols",
            limitations=_LIMITATIONS,
            reason=screen.reason,
        )

    volume = profile.dose_mg / math.exp(screen.intercept)
    clearance = screen.kel * volume
    predicted = tuple(
        (profile.dose_mg / volume) * math.exp(-(clearance / volume) * point.time_h) for point in profile.points
    )
    relative = [
        abs(pred - point.concentration_mg_l) / point.concentration_mg_l
        for pred, point in zip(predicted, profile.points, strict=True)
    ]
    max_relative = max(relative)
    status = DiagnosticStatus.PASS if max_relative <= 0.01 else DiagnosticStatus.WARNING
    r_squared = 0.0 if screen.r_squared is None else screen.r_squared
    return _succeeded(
        profile,
        method="log_linear_ols",
        limitations=_LIMITATIONS,
        status=status,
        clearance=clearance,
        volume=volume,
        predicted=predicted,
        warnings=[] if status is DiagnosticStatus.PASS else ["A relative residual exceeds 1%."],
        metrics=[
            Metric(name="log_r_squared", value=r_squared),
            Metric(name="max_relative_residual", value=max_relative),
            Metric(name="n_points", value=float(len(profile.points))),
        ],
    )


def screen_log_linear(points: tuple[QuantifiedPoint, ...] | list[QuantifiedPoint]) -> LogLinearScreen:
    """Screen quantified points for a declining log-linear IV bolus.

    Fewer than three points, a non-positive concentration, or times that do
    not vary are not identifiable. A non-negative slope is identifiable as a
    regression but is not a declining concentration, so clearance is withheld.
    """

    series = tuple(points)
    if len(series) < 3:
        return LogLinearScreen(False, False, None, None, None, None, "Fewer than three quantified points.")
    if any(point.concentration_mg_l <= 0 for point in series):
        return LogLinearScreen(
            False,
            False,
            None,
            None,
            None,
            None,
            "A quantified concentration is not positive.",
        )
    xs = [point.time_h for point in series]
    ys = [math.log(point.concentration_mg_l) for point in series]
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0:
        return LogLinearScreen(
            False,
            False,
            None,
            None,
            None,
            None,
            "Quantified sample times do not vary, so CL and V are not identifiable.",
        )
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    intercept = ybar - slope * xbar
    if not math.isfinite(slope) or not math.isfinite(intercept):
        return LogLinearScreen(
            False,
            False,
            None,
            None,
            None,
            None,
            "The log-linear regression is not finite, so CL and V are not identifiable.",
        )
    sst = sum((y - ybar) ** 2 for y in ys)
    sse = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys, strict=True))
    r_squared = 1.0 if sst <= 1e-18 else 1.0 - max(0.0, sse) / sst
    kel = -slope
    if kel <= 0 or math.exp(intercept) <= 0:
        return LogLinearScreen(
            True,
            False,
            None,
            None,
            slope,
            r_squared,
            "The log-linear slope is not a declining concentration.",
        )
    return LogLinearScreen(True, True, kel, intercept, slope, r_squared, "")


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
) -> IvBolusFit:
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
        model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_IV_BOLUS.context_of_use,
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
            route=Route.IV_BOLUS,
        ),
        ParameterEstimate(
            fit_id=fit_id,
            parameter="V",
            semantic=ParameterSemantic.V,
            value=volume,
            unit="L",
            estimation_method=method,
            route=Route.IV_BOLUS,
        ),
    )
    return _pack(profile, fit, estimates, diagnostic, predicted)


def _failed(
    profile: ConcentrationProfile,
    *,
    method: str,
    limitations: tuple[str, ...],
    reason: str,
) -> IvBolusFit:
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
        model_id=ONE_COMPARTMENT_IV_BOLUS.model_id,
        dataset_id=_dataset_id(profile),
        software_version=__version__,
        estimation_method=method,
        context_of_use=ONE_COMPARTMENT_IV_BOLUS.context_of_use,
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
) -> IvBolusFit:
    return IvBolusFit(
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
