"""Direct-effect Emax with a fixed baseline and Hill coefficient 1.

The estimating equation is the linear rearrangement of

    E = E0 + Emax * C / (EC50 + C)

on quantified pairs. It does not estimate a Hill coefficient, a lag, or an
indirect response. Predictions are labeled predicted and are not effect
observations.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, PositiveFloat, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.modeling import DiagnosticRecord, DiagnosticStatus, Metric, ModelFit
from pharmacoweave.core.units import CONC_AMOUNT, CONC_MASS, convert, parse_unit
from pharmacoweave.errors import AnalysisError, ContractViolation, UnitError
from pharmacoweave.pd_models.records import EffectDirection
from pharmacoweave.pd_models.specs import DIRECT_EMAX
from pharmacoweave.version import __version__

_METHOD = "direct_emax_linearized_v1"
_LIMITATIONS = (
    RESEARCH_USE_LIMITATION,
    "Direct Emax with a fixed baseline and hill fixed at 1. No standard error is reported.",
    "Each effect is paired with an explicit concentration. AUC, Cmax, lag, and indirect response are not used.",
    "A failed diagnostic does not report EC50.",
)


class EffectParameter(Record):
    """One direct-effect parameter. EC50 is a concentration, not a clearance."""

    fit_id: Identifier
    name: str = Field(min_length=1)
    value: FiniteFloat
    unit: str
    fixed: bool
    estimation_method: str = Field(min_length=1)

    @model_validator(mode="after")
    def unit_matches_name(self) -> "EffectParameter":
        parsed = parse_unit(self.unit)
        if self.name == "EC50":
            if self.value <= 0:
                raise ValueError("EC50 must be positive.")
            if parsed.dimensions not in {CONC_MASS, CONC_AMOUNT}:
                raise ValueError("EC50 must use a concentration unit.")
        if self.name == "hill":
            if self.value <= 0:
                raise ValueError("hill must be positive.")
            if self.unit != "1":
                raise ValueError("hill uses the dimensionless unit 1.")
        return self


class DirectEffectPair(Record):
    """One observed effect and the concentration it is paired with."""

    effect_id: Identifier
    exposure_ref: Identifier
    concentration: PositiveFloat
    concentration_unit: str
    effect: FiniteFloat
    effect_unit: str

    @model_validator(mode="after")
    def units(self) -> "DirectEffectPair":
        if parse_unit(self.concentration_unit).dimensions not in {CONC_MASS, CONC_AMOUNT}:
            raise ValueError("concentration_unit must be a concentration.")
        parse_unit(self.effect_unit)
        return self


@dataclass(frozen=True)
class PredictedEffect:
    """A model effect. ``evidence_grade`` stays predicted."""

    effect_id: str
    exposure_ref: str
    value: float
    unit: str
    derived_from: tuple[str, ...]
    evidence_grade: str = "predicted"

    def __post_init__(self) -> None:
        if self.evidence_grade != "predicted":
            raise ValueError("A model effect is predicted and cannot be labeled observed.")


@dataclass(frozen=True)
class DirectEmaxFit:
    """One direct Emax fit. Failed fits have empty estimates and predictions."""

    fit: ModelFit
    estimates: tuple[EffectParameter, ...]
    diagnostic: DiagnosticRecord
    predicted: tuple[PredictedEffect, ...]
    effect_ids: tuple[str, ...]
    exposure_refs: tuple[str, ...]
    reason: str


def emax_effect(*, concentration: float, e0: float, emax: float, ec50: float, hill: float = 1.0) -> float:
    """Return E0 + Emax * C^hill / (EC50^hill + C^hill).

    This is a forward evaluation. ``fit_direct_emax`` estimates only the
    hill=1 case. ``estimate_hill`` refuses to estimate hill.
    """

    if concentration < 0 or ec50 <= 0 or hill <= 0 or not math.isfinite(e0) or not math.isfinite(emax):
        raise ValueError("concentration must be >= 0, EC50 and hill must be > 0, and E0 and Emax must be finite.")
    if concentration == 0:
        return e0
    powered = concentration**hill
    return e0 + emax * powered / (ec50**hill + powered)


def estimate_hill() -> None:
    """Refuse to estimate a Hill coefficient."""

    raise AnalysisError(
        "Sigmoid Emax hill-coefficient estimation is not implemented. "
        "fit_direct_emax fixes hill at 1 and does not estimate a delay or an indirect response."
    )


def effect_observation_from_prediction(predicted: PredictedEffect) -> None:
    """Refuse to store a prediction as an observed effect."""

    raise ContractViolation(
        f"{predicted.effect_id} is a model prediction and cannot be stored as an observed effect."
    )


def fit_direct_emax(
    pairs: Sequence[DirectEffectPair],
    *,
    e0: float,
    e0_unit: str,
    direction: EffectDirection = EffectDirection.UNSPECIFIED,
    fit_id: str = "direct-emax",
) -> DirectEmaxFit:
    """Estimate Emax and EC50. Non-identifiable curves fail without an EC50."""

    if not math.isfinite(e0):
        raise AnalysisError("e0 must be finite.")
    parse_unit(e0_unit)
    rows = tuple(pairs)
    effect_ids = tuple(row.effect_id for row in rows)
    exposure_refs = tuple(row.exposure_ref for row in rows)
    if len(effect_ids) != len(set(effect_ids)):
        raise AnalysisError("effect_id values must be unique inside one direct Emax fit.")
    if len(exposure_refs) != len(set(exposure_refs)):
        raise AnalysisError("exposure_ref values must be unique inside one direct Emax fit.")
    if len(rows) < 3:
        return _failed(fit_id, effect_ids, exposure_refs, "Fewer than three quantified concentration-effect pairs.")

    concentration_unit = rows[0].concentration_unit
    try:
        concentrations = [convert(row.concentration, row.concentration_unit, concentration_unit) for row in rows]
        effects = [convert(row.effect, row.effect_unit, e0_unit) for row in rows]
    except UnitError as exc:
        raise AnalysisError(str(exc)) from exc

    scale = max(1.0, abs(e0), max(abs(value) for value in effects))
    if any(abs(value - e0) <= 1e-12 * scale for value in effects):
        return _failed(
            fit_id,
            effect_ids,
            exposure_refs,
            "An effect equals the fixed baseline, so the linear Emax rearrangement is undefined.",
        )
    deltas = [value - e0 for value in effects]
    if any(delta > 0 for delta in deltas) and any(delta < 0 for delta in deltas):
        return _failed(
            fit_id,
            effect_ids,
            exposure_refs,
            "Effects move both above and below the fixed baseline, so one Emax is not identifiable.",
        )
    span = max(deltas) - min(deltas)
    if span <= 1e-12 * max(1.0, max(abs(delta) for delta in deltas)):
        return _failed(
            fit_id,
            effect_ids,
            exposure_refs,
            "Effects do not vary with concentration, so Emax and EC50 are not identifiable.",
        )
    if max(concentrations) - min(concentrations) <= 0:
        return _failed(fit_id, effect_ids, exposure_refs, "Concentrations do not vary, so EC50 is not identifiable.")

    slope_intercept = _line([1.0 / value for value in concentrations], [1.0 / delta for delta in deltas])
    if slope_intercept is None:
        return _failed(fit_id, effect_ids, exposure_refs, "Concentrations do not vary, so EC50 is not identifiable.")
    slope, intercept = slope_intercept
    if not math.isfinite(slope) or not math.isfinite(intercept) or intercept == 0:
        return _failed(fit_id, effect_ids, exposure_refs, "Emax is not identifiable.")
    emax = 1.0 / intercept
    ec50 = slope / intercept
    if not math.isfinite(emax) or not math.isfinite(ec50) or ec50 <= 0:
        return _failed(fit_id, effect_ids, exposure_refs, "EC50 is not positive, so the direct Emax fit is rejected.")
    if direction is EffectDirection.INCREASE and emax <= 0:
        return _failed(
            fit_id,
            effect_ids,
            exposure_refs,
            "The declared direction is increase, but Emax is not positive. EC50 was not reported.",
        )
    if direction is EffectDirection.DECREASE and emax >= 0:
        return _failed(
            fit_id,
            effect_ids,
            exposure_refs,
            "The declared direction is decrease, but Emax is not negative. EC50 was not reported.",
        )

    predicted_values = [
        emax_effect(concentration=value, e0=e0, emax=emax, ec50=ec50, hill=1.0) for value in concentrations
    ]
    residual = max(abs(pred - observed) for pred, observed in zip(predicted_values, effects, strict=True))
    status = DiagnosticStatus.PASS if residual <= 1e-6 * scale else DiagnosticStatus.WARNING
    warnings = [] if status is DiagnosticStatus.PASS else ["An effect residual exceeds the noise-free tolerance."]
    estimates = (
        EffectParameter(
            fit_id=fit_id,
            name="E0",
            value=e0,
            unit=e0_unit,
            fixed=True,
            estimation_method="fixed",
        ),
        EffectParameter(
            fit_id=fit_id,
            name="Emax",
            value=emax,
            unit=e0_unit,
            fixed=False,
            estimation_method=_METHOD,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="EC50",
            value=ec50,
            unit=concentration_unit,
            fixed=False,
            estimation_method=_METHOD,
        ),
        EffectParameter(
            fit_id=fit_id,
            name="hill",
            value=1.0,
            unit="1",
            fixed=True,
            estimation_method="fixed",
        ),
    )
    predicted = tuple(
        PredictedEffect(
            effect_id=row.effect_id,
            exposure_ref=row.exposure_ref,
            value=value,
            unit=e0_unit,
            derived_from=(f"fit:{fit_id}", f"exposure:{row.exposure_ref}"),
        )
        for row, value in zip(rows, predicted_values, strict=True)
    )
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=status,
        limitations=list(_LIMITATIONS),
        warnings=warnings,
        metrics=[
            Metric(name="n_points", value=float(len(rows))),
            Metric(name="max_abs_residual", value=residual, unit=e0_unit),
            Metric(name="Emax", value=emax, unit=e0_unit),
            Metric(name="EC50", value=ec50, unit=concentration_unit),
        ],
    )
    fit = _model_fit(fit_id, status, converged=True, limitations=list(_LIMITATIONS))
    return DirectEmaxFit(fit, estimates, diagnostic, predicted, effect_ids, exposure_refs, "")


def _line(xs: list[float], ys: list[float]) -> tuple[float, float] | None:
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0:
        return None
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=True)) / sxx
    return slope, ybar - slope * xbar


def _failed(
    fit_id: str,
    effect_ids: tuple[str, ...],
    exposure_refs: tuple[str, ...],
    reason: str,
) -> DirectEmaxFit:
    stated = [*_LIMITATIONS, reason]
    diagnostic = DiagnosticRecord(
        fit_id=fit_id,
        status=DiagnosticStatus.FAIL,
        limitations=stated,
        warnings=["Emax and EC50 were not estimated."],
    )
    fit = _model_fit(fit_id, DiagnosticStatus.FAIL, converged=False, limitations=stated)
    return DirectEmaxFit(fit, (), diagnostic, (), effect_ids, exposure_refs, reason)


def _model_fit(fit_id: str, status: DiagnosticStatus, *, converged: bool, limitations: list[str]) -> ModelFit:
    return ModelFit(
        fit_id=fit_id,
        model_id=DIRECT_EMAX.model_id,
        dataset_id=fit_id,
        software_version=__version__,
        estimation_method=_METHOD,
        context_of_use=DIRECT_EMAX.context_of_use,
        diagnostic_status=status,
        converged=converged,
        limitations=limitations,
    )
