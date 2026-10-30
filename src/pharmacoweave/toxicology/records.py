"""Research-use toxicology endpoint records.

A record can cite a dose event or an exposure id. It cannot be turned into
a clinical dose, a safety margin, or a treatment recommendation.
"""

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, Record
from pharmacoweave.core.evidence import TOX_RESEARCH_LIMITATION
from pharmacoweave.core.units import (
    AMOUNT,
    AMOUNT_PER_BODY_MASS,
    BODY_MASS,
    CONC_AMOUNT,
    CONC_MASS,
    DIMENSIONLESS,
    DOSE_PER_BODY_MASS,
    DRUG_MASS,
    RATE,
    parse_unit,
)
from pharmacoweave.errors import ContractViolation


class ToxKind(StrEnum):
    BODY_WEIGHT_CHANGE = "body_weight_change"
    IN_LIFE_OBSERVATION = "in_life_observation"
    SERUM_CHEMISTRY = "serum_chemistry"
    HISTOPATHOLOGY = "histopathology"
    MORTALITY = "mortality"


class ToxSeverity(StrEnum):
    NOT_GRADED = "not_graded"
    MINIMAL = "minimal"
    MILD = "mild"
    MODERATE = "moderate"
    MARKED = "marked"


_DOSE_DIMENSIONS = {DRUG_MASS, AMOUNT, DOSE_PER_BODY_MASS, AMOUNT_PER_BODY_MASS}
_PROPORTION = {ToxKind.IN_LIFE_OBSERVATION, ToxKind.HISTOPATHOLOGY, ToxKind.MORTALITY}


class ToxEndpoint(Record):
    """One observed toxicology finding. ``evidence_grade`` stays observed."""

    endpoint_id: Identifier
    compound_id: Identifier
    kind: ToxKind
    finding: str = Field(min_length=1)
    severity: ToxSeverity
    limitations: list[str] = Field(min_length=1)
    subject_id: Identifier | None = None
    dose_event_id: Identifier | None = None
    exposure_ref: str | None = None
    numeric_value: FiniteFloat | None = None
    unit: str | None = None
    incidence_affected: int | None = Field(default=None, ge=0)
    incidence_examined: int | None = Field(default=None, ge=0)
    evidence_grade: Literal["observed"] = "observed"
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def links_and_units(self) -> "ToxEndpoint":
        if TOX_RESEARCH_LIMITATION not in self.limitations:
            raise ValueError(f"limitations must include {TOX_RESEARCH_LIMITATION!r}.")
        if self.dose_event_id is None and not self.exposure_ref:
            raise ValueError("A toxicology endpoint must reference a dose_event_id or an exposure_ref.")
        if self.exposure_ref is not None and not self.exposure_ref.strip():
            raise ValueError("exposure_ref cannot be blank.")
        if (self.numeric_value is None) != (self.unit is None):
            raise ValueError("numeric_value and unit are recorded together.")
        if (self.incidence_affected is None) != (self.incidence_examined is None):
            raise ValueError("incidence_affected and incidence_examined are recorded together.")
        if self.incidence_examined is not None:
            if self.incidence_examined < 1:
                raise ValueError("incidence_examined must be at least 1.")
            if self.incidence_affected is not None and self.incidence_affected > self.incidence_examined:
                raise ValueError("incidence_affected cannot exceed incidence_examined.")
        if self.kind is ToxKind.MORTALITY and self.incidence_affected is None:
            raise ValueError("Mortality requires an incidence.")
        if self.kind is ToxKind.HISTOPATHOLOGY and self.severity is ToxSeverity.NOT_GRADED:
            raise ValueError("Histopathology requires a severity grade.")
        if self.numeric_value is not None and self.unit is not None:
            _check_unit(self.kind, self.unit)
        return self


def dosing_recommendation(record: ToxEndpoint) -> None:
    """Refuse to turn a toxicology observation into a clinical dose."""

    raise ContractViolation(
        f"{record.endpoint_id} is a research-use toxicology observation. "
        "It does not produce a clinical dose, a safety margin, or a treatment recommendation."
    )


def _check_unit(kind: ToxKind, unit: str) -> None:
    parsed = parse_unit(unit)
    if parsed.dimensions in _DOSE_DIMENSIONS:
        raise ValueError("A toxicology endpoint cannot be recorded as an administered dose.")
    if kind is ToxKind.BODY_WEIGHT_CHANGE:
        proportion = parsed.dimensions == DIMENSIONLESS and parsed.factor in {1.0, 0.01}
        if parsed.dimensions != BODY_MASS and not proportion:
            raise ValueError("body_weight_change uses body mass, percent, or a dimensionless fraction.")
        return
    if kind is ToxKind.SERUM_CHEMISTRY:
        allowed = parsed.dimensions in {CONC_MASS, CONC_AMOUNT, RATE} or (
            parsed.dimensions == DIMENSIONLESS and parsed.factor in {1.0, 0.01}
        )
        if not allowed:
            raise ValueError("serum_chemistry uses a concentration, rate, percent, or dimensionless fraction.")
        return
    if kind in _PROPORTION:
        if not (parsed.dimensions == DIMENSIONLESS and parsed.factor in {1.0, 0.01}):
            raise ValueError(f"{kind.value} numeric results use percent or a dimensionless fraction.")
