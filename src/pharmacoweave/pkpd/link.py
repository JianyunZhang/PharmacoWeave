"""Link an observed direct effect to an explicit concentration.

AUC and Cmax remain exposure summaries. They are not rewritten as the
concentration driver of a direct Emax model. BLQ effects are omitted and are
not replaced with zero.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

from pydantic import model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, ObservationStatus, Record
from pharmacoweave.core.units import AUC_AMOUNT, AUC_MASS, CONC_AMOUNT, CONC_MASS, parse_unit
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pd_models.emax import DirectEffectPair, DirectEmaxFit, fit_direct_emax
from pharmacoweave.pd_models.records import EffectDirection, EffectObservation, EffectTimeScale


class ExposureKind(StrEnum):
    CONCENTRATION = "concentration"
    AUC = "auc"
    CMAX = "cmax"


class ExposureValue(Record):
    """One exposure number. Kind and unit must agree."""

    source_id: Identifier
    kind: ExposureKind
    value: FiniteFloat
    unit: str

    @model_validator(mode="after")
    def kind_and_unit(self) -> "ExposureValue":
        parsed = parse_unit(self.unit)
        if self.kind in {ExposureKind.CONCENTRATION, ExposureKind.CMAX}:
            if parsed.dimensions not in {CONC_MASS, CONC_AMOUNT}:
                raise ValueError("concentration and cmax require a concentration unit.")
            if self.value <= 0:
                raise ValueError("concentration and cmax must be positive.")
        elif parsed.dimensions not in {AUC_MASS, AUC_AMOUNT}:
            raise ValueError("auc requires an AUC unit.")
        elif self.value <= 0:
            raise ValueError("auc must be positive.")
        return self


@dataclass(frozen=True)
class OmittedEffect:
    """An effect that was not sent to the direct Emax fit."""

    effect_id: str
    reason: str


@dataclass(frozen=True)
class LinkedDirectEffects:
    """Quantified direct pairs and the effects that were left out."""

    pairs: tuple[DirectEffectPair, ...]
    omitted: tuple[OmittedEffect, ...]
    direction: EffectDirection


def fit_effect_lag() -> None:
    """Refuse to estimate a delay between concentration and effect."""

    raise AnalysisError(
        "Effect lag is not estimated. A direct effect records time_since_exposure but does not fit a delay parameter."
    )


def fit_indirect_response() -> None:
    """Refuse to fit an indirect-response model."""

    raise AnalysisError(
        "Indirect-response PK/PD is not implemented. This version links a direct effect to an explicit concentration."
    )


def link_direct_effects(
    effects: Sequence[EffectObservation],
    exposures: Sequence[ExposureValue],
) -> LinkedDirectEffects:
    """Pair quantified direct effects with concentrations.

    Non-quantified rows and non-direct time scales are omitted. A quantified
    direct effect that points at AUC or Cmax is an error, not a silent
    substitution.
    """

    by_id: dict[str, ExposureValue] = {}
    for exposure in exposures:
        if exposure.source_id in by_id:
            raise AnalysisError(f"Duplicate exposure source_id {exposure.source_id}.")
        by_id[exposure.source_id] = exposure

    seen: set[str] = set()
    pairs: list[DirectEffectPair] = []
    omitted: list[OmittedEffect] = []
    directions: list[EffectDirection] = []
    for effect in effects:
        if effect.effect_id in seen:
            raise AnalysisError(f"Duplicate effect_id {effect.effect_id}.")
        seen.add(effect.effect_id)
        if effect.status is not ObservationStatus.QUANTIFIED or effect.numeric_value is None or effect.unit is None:
            omitted.append(
                OmittedEffect(
                    effect.effect_id,
                    f"status is {effect.status.value}; the effect was not replaced with zero.",
                )
            )
            continue
        if effect.time_scale is not EffectTimeScale.DIRECT:
            omitted.append(
                OmittedEffect(
                    effect.effect_id,
                    f"time_scale is {effect.time_scale.value}; only direct effects are linked.",
                )
            )
            continue
        matched = by_id.get(effect.exposure_ref)
        if matched is None:
            raise AnalysisError(f"No exposure is recorded for {effect.effect_id} ({effect.exposure_ref}).")
        if matched.kind == ExposureKind.AUC:
            raise AnalysisError("AUC is an exposure summary, not a direct-effect concentration.")
        if matched.kind == ExposureKind.CMAX:
            raise AnalysisError("Cmax is an exposure summary, not a direct-effect concentration.")
        pairs.append(
            DirectEffectPair(
                effect_id=effect.effect_id,
                exposure_ref=effect.exposure_ref,
                concentration=matched.value,
                concentration_unit=matched.unit,
                effect=effect.numeric_value,
                effect_unit=effect.unit,
            )
        )
        directions.append(effect.direction)

    if not directions:
        direction = EffectDirection.UNSPECIFIED
    elif any(item is not directions[0] for item in directions):
        raise AnalysisError("Direct effects in one link do not share a direction.")
    else:
        direction = directions[0]
    return LinkedDirectEffects(tuple(pairs), tuple(omitted), direction)


def fit_linked_direct_emax(
    effects: Sequence[EffectObservation],
    exposures: Sequence[ExposureValue],
    *,
    e0: float,
    e0_unit: str,
    fit_id: str = "direct-emax",
) -> tuple[LinkedDirectEffects, DirectEmaxFit]:
    """Link effects to concentrations and fit direct Emax."""

    linked = link_direct_effects(effects, exposures)
    fit = fit_direct_emax(
        linked.pairs,
        e0=e0,
        e0_unit=e0_unit,
        direction=linked.direction,
        fit_id=fit_id,
    )
    return linked, fit
