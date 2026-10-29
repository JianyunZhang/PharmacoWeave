"""Experimental entities for a pharmacology study.

Single-record rules live on the models. Foreign keys and clock checks live
in ``pharmacoweave.core.validation``.
"""

import math
from datetime import UTC, datetime
from enum import StrEnum
from typing import Annotated, Literal

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, model_validator

from pharmacoweave.core.units import (
    AMOUNT,
    BODY_MASS,
    CONC_AMOUNT,
    CONC_MASS,
    DRUG_MASS,
    MASS_PER_AMOUNT,
    TIME,
    parse_unit,
)

_ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$"


def _finite(value: float) -> float:
    if not math.isfinite(value):
        raise ValueError("Value must be finite.")
    return value


def _aware_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError("Absolute time must include a timezone offset.")
    return value.astimezone(UTC)


FiniteFloat = Annotated[float, AfterValidator(_finite)]
PositiveFloat = Annotated[float, Field(gt=0), AfterValidator(_finite)]
Identifier = Annotated[str, Field(pattern=_ID_PATTERN)]
UtcDateTime = Annotated[datetime, AfterValidator(_aware_utc)]


class Record(BaseModel):
    """Base for immutable records that reject unknown fields."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class StudyType(StrEnum):
    PK = "pk"
    PD = "pd"
    PKPD = "pkpd"
    TK = "tk"
    DMPK = "dmpk"
    TOXICOLOGY = "toxicology"
    OTHER = "other"


class DataClass(StrEnum):
    SYNTHETIC = "synthetic"
    PRECLINICAL_ANIMAL = "preclinical_animal"
    CLINICAL = "clinical"
    LITERATURE = "literature"
    OTHER = "other"


class Sex(StrEnum):
    FEMALE = "female"
    MALE = "male"
    UNKNOWN = "unknown"
    NOT_APPLICABLE = "not_applicable"


class AnalyteRole(StrEnum):
    PARENT = "parent"
    METABOLITE = "metabolite"
    BIOMARKER = "biomarker"
    OTHER = "other"


class AnalyteForm(StrEnum):
    TOTAL = "total"
    FREE = "free"
    PARENT = "parent"
    METABOLITE = "metabolite"
    OCCUPANCY = "occupancy"
    UNSPECIFIED = "unspecified"


class Route(StrEnum):
    IV_BOLUS = "iv_bolus"
    IV_INFUSION = "iv_infusion"
    ORAL = "oral"
    SUBCUTANEOUS = "subcutaneous"
    INTRAMUSCULAR = "intramuscular"
    OTHER = "other"


class ObservationType(StrEnum):
    CONCENTRATION = "concentration"
    EFFECT = "effect"
    TOXICITY = "toxicity"


class ObservationStatus(StrEnum):
    QUANTIFIED = "quantified"
    BLQ = "BLQ"
    ABOVE_ULOQ = "above_ULOQ"
    MISSING = "missing"
    INVALID = "invalid"


class AssayQuantitation(StrEnum):
    CONCENTRATION = "concentration"
    OTHER = "other"


class QcStatus(StrEnum):
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    NOT_REVIEWED = "not_reviewed"


class Quantity(Record):
    """A finite number with an explicit unit."""

    value: FiniteFloat
    unit: str

    @model_validator(mode="after")
    def unit_parses(self) -> "Quantity":
        parse_unit(self.unit)
        return self


class Timebase(Record):
    """A named clock. Samples declare which clock their times use."""

    timebase_id: Identifier
    description: str = Field(min_length=1)


class StudyMetadata(Record):
    """Design context for one study. This is not an analysis result."""

    study_id: Identifier
    title: str = Field(min_length=1)
    study_type: StudyType
    data_class: DataClass
    species: str = Field(min_length=1)
    design_version: str = Field(min_length=1)
    timebases: list[Timebase] = Field(min_length=1)
    protocol_id: str | None = None
    ncbi_taxon_id: int | None = Field(default=None, gt=0)
    description: str | None = None
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def timebases_are_unique(self) -> "StudyMetadata":
        ids = [item.timebase_id for item in self.timebases]
        if len(ids) != len(set(ids)):
            raise ValueError("timebase_id values must be unique inside a study.")
        return self


class Subject(Record):
    """One experimental unit. Technical replicates are not subjects."""

    subject_id: Identifier
    study_id: Identifier
    species: str = Field(min_length=1)
    sex: Sex = Sex.UNKNOWN
    group_id: str | None = None
    body_weight: Quantity | None = None

    @model_validator(mode="after")
    def body_weight_dimension(self) -> "Subject":
        if self.body_weight is not None:
            dimensions = parse_unit(self.body_weight.unit).dimensions
            if dimensions != BODY_MASS:
                raise ValueError("body_weight must use a body-mass unit such as kg.")
            if self.body_weight.value <= 0:
                raise ValueError("body_weight must be positive.")
        return self


class Compound(Record):
    """A test article. A TheraWeave id is a cross-reference, not a measurement."""

    compound_id: Identifier
    name: str = Field(min_length=1)
    molecular_weight: PositiveFloat | None = None
    molecular_weight_unit: str | None = None
    structure_ref: str | None = None
    formulation_ref: str | None = None
    theraweave_candidate_id: str | None = None

    @model_validator(mode="after")
    def molecular_weight_pair(self) -> "Compound":
        present = (self.molecular_weight is None, self.molecular_weight_unit is None)
        if present == (True, False) or present == (False, True):
            raise ValueError("molecular_weight and molecular_weight_unit are recorded together.")
        if self.molecular_weight_unit is not None:
            dimensions = parse_unit(self.molecular_weight_unit).dimensions
            if dimensions != MASS_PER_AMOUNT:
                raise ValueError("molecular_weight_unit must be drug mass per amount, such as g/mol.")
        return self


class Analyte(Record):
    """What was measured. Parent, metabolite, and biomarker stay distinct."""

    analyte_id: Identifier
    name: str = Field(min_length=1)
    analyte_role: AnalyteRole
    matrix: str = Field(min_length=1)
    compound_id: Identifier | None = None

    @model_validator(mode="after")
    def compound_required_for_drug_analyte(self) -> "Analyte":
        if self.analyte_role in {AnalyteRole.PARENT, AnalyteRole.METABOLITE} and self.compound_id is None:
            raise ValueError("Parent and metabolite analytes require compound_id.")
        return self


class Assay(Record):
    """Bioanalytical method limits. LLOQ is required for a concentration assay."""

    assay_id: Identifier
    method: str = Field(min_length=1)
    analyte_id: Identifier
    matrix: str = Field(min_length=1)
    quantitation: AssayQuantitation = AssayQuantitation.CONCENTRATION
    qc_status: QcStatus = QcStatus.NOT_REVIEWED
    calibration_id: str | None = None
    lloq: PositiveFloat | None = None
    uloq: PositiveFloat | None = None
    limit_unit: str | None = None

    @model_validator(mode="after")
    def concentration_limits(self) -> "Assay":
        if self.quantitation != AssayQuantitation.CONCENTRATION:
            return self
        if self.lloq is None or self.limit_unit is None:
            raise ValueError("A concentration assay requires lloq and limit_unit.")
        dimensions = parse_unit(self.limit_unit).dimensions
        if dimensions not in {CONC_MASS, CONC_AMOUNT}:
            raise ValueError("limit_unit must be a mass or amount concentration.")
        if self.uloq is not None and self.uloq <= self.lloq:
            raise ValueError("uloq must be greater than lloq.")
        return self


class DoseEvent(Record):
    """An administered dose. The primary amount is the actual dose, not the plan."""

    dose_event_id: Identifier
    subject_id: Identifier
    compound_id: Identifier
    absolute_time: UtcDateTime
    dose_value: PositiveFloat
    dose_unit: str
    route: Route
    formulation: str | None = None
    infusion_duration: PositiveFloat | None = None
    infusion_duration_unit: str | None = None
    body_weight_normalized_value: PositiveFloat | None = None
    body_weight_normalized_unit: str | None = None

    @model_validator(mode="after")
    def dose_units_and_route(self) -> "DoseEvent":
        dimensions = parse_unit(self.dose_unit).dimensions
        if dimensions not in {DRUG_MASS, AMOUNT}:
            raise ValueError(
                "dose_value must be an administered amount, such as mg or nmol. "
                "Record a per-body-weight dose in the normalized fields."
            )
        normalized = (self.body_weight_normalized_value is None, self.body_weight_normalized_unit is None)
        if normalized in {(True, False), (False, True)}:
            raise ValueError("Normalized dose value and unit are recorded together.")
        if self.body_weight_normalized_unit is not None:
            normalized_dimensions = parse_unit(self.body_weight_normalized_unit).dimensions
            allowed = {
                DRUG_MASS: (0, 1, -1, 0, 0),
                AMOUNT: (0, 0, -1, 0, 1),
            }[dimensions]
            if normalized_dimensions != allowed:
                raise ValueError("Normalized dose unit does not match the administered amount.")
        infusion = (self.infusion_duration is None, self.infusion_duration_unit is None)
        if self.route == Route.IV_INFUSION:
            if infusion != (False, False):
                raise ValueError("IV infusion requires infusion_duration and infusion_duration_unit.")
            if parse_unit(self.infusion_duration_unit or "").dimensions != TIME:
                raise ValueError("infusion_duration_unit must be a time unit.")
        elif infusion != (True, True):
            raise ValueError("infusion_duration is recorded only for route iv_infusion.")
        return self


class Biosample(Record):
    """A collected sample tied to one subject, one clock, and one dose event."""

    sample_id: Identifier
    subject_id: Identifier
    collection_time_absolute: UtcDateTime
    time_since_dose: FiniteFloat
    time_since_dose_unit: str
    dose_event_id: Identifier
    matrix: str = Field(min_length=1)
    timebase_id: Identifier
    is_predose: bool = False
    collection_site: str | None = None

    @model_validator(mode="after")
    def relative_time(self) -> "Biosample":
        if parse_unit(self.time_since_dose_unit).dimensions != TIME:
            raise ValueError("time_since_dose_unit must be a time unit.")
        if self.time_since_dose < 0 and not self.is_predose:
            raise ValueError("A negative time_since_dose is predose and requires is_predose=true.")
        if self.is_predose and self.time_since_dose >= 0:
            raise ValueError("is_predose requires a negative time_since_dose.")
        return self


class Observation(Record):
    """One measured result. BLQ, missing, and invalid rows do not carry a number."""

    observation_id: Identifier
    sample_id: Identifier
    observation_type: ObservationType
    status: ObservationStatus
    analyte_form: AnalyteForm = AnalyteForm.UNSPECIFIED
    analyte_id: Identifier | None = None
    raw_value: str | None = None
    reported_value: str | None = None
    numeric_value: FiniteFloat | None = None
    unit: str | None = None
    assay_id: Identifier | None = None
    qc_flags: list[str] = Field(default_factory=list)
    exclusion_reason: str | None = None
    replicate_id: str | None = None
    assay_batch_id: str | None = None

    @model_validator(mode="after")
    def status_and_value(self) -> "Observation":
        for flag in self.qc_flags:
            if not flag or "|" in flag:
                raise ValueError("qc_flags entries must be non-empty and cannot contain '|'.")
        if self.observation_type == ObservationType.CONCENTRATION and self.analyte_id is None:
            raise ValueError("A concentration observation requires analyte_id.")
        if self.observation_type == ObservationType.CONCENTRATION and self.assay_id is None:
            raise ValueError("A concentration observation requires assay_id.")
        if self.unit is not None:
            dimensions = parse_unit(self.unit).dimensions
            if self.observation_type == ObservationType.CONCENTRATION and dimensions not in {
                CONC_MASS,
                CONC_AMOUNT,
            }:
                raise ValueError("A concentration unit must be mass/volume or amount/volume.")
        if self.status == ObservationStatus.QUANTIFIED:
            if self.numeric_value is None or self.unit is None:
                raise ValueError("A quantified observation requires numeric_value and unit.")
            if self.observation_type == ObservationType.CONCENTRATION and self.numeric_value <= 0:
                raise ValueError("A quantified concentration must be positive. BLQ is not recorded as zero.")
        elif self.numeric_value is not None:
            raise ValueError(
                "numeric_value is recorded only when status is quantified. "
                "BLQ, above_ULOQ, missing, and invalid stay non-numeric."
            )
        if self.status == ObservationStatus.INVALID and not self.exclusion_reason:
            raise ValueError("An invalid observation requires exclusion_reason.")
        return self
