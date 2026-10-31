"""In vitro ADME records.

Permeability, microsomal clearance, hepatocyte clearance, plasma and
microsomal unbound fractions, and blood-to-plasma ratios are stored as
reported in vitro results. A plasma clearance unit is rejected.
``scale_to_in_vivo_clearance`` still refuses to turn one record into an in vivo
CL; the exploratory IVIVE in ``pharmacoweave.dmpk.ivive`` requires every
physiological assumption to be passed explicitly and returns a predicted
record, never an observed or estimated CL.
"""

from collections.abc import Sequence
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.units import CLEARANCE, parse_unit
from pharmacoweave.dmpk.physiology import IviveClearancePrediction, PhysiologicalScalingParams
from pharmacoweave.errors import AnalysisError, UnitError
from pharmacoweave.version import SCHEMA_VERSION

PERMEABILITY_TO_CM_S = {"cm/s": 1.0, "um/s": 1e-4, "nm/s": 1e-7}
MICROSOMAL_TO_UL_MIN_MG = {"uL/min/mg_protein": 1.0, "mL/min/mg_protein": 1000.0}
HEPATOCYTE_TO_UL_MIN_MCELL = {"uL/min/million_cells": 1.0, "mL/min/million_cells": 1000.0}
UNBOUND_FRACTION_TO_FRACTION = {"fraction": 1.0, "1": 1.0, "percent": 0.01, "pct": 0.01}
BLOOD_TO_PLASMA_TO_RATIO = {"ratio": 1.0, "1": 1.0, "fold": 1.0}


class AdmeAssayKind(StrEnum):
    PERMEABILITY = "permeability"
    MICROSOMAL_CLEARANCE = "microsomal_clearance"
    HEPATOCYTE_CLEARANCE = "hepatocyte_clearance"
    FU_PLASMA = "fu_plasma"
    FU_MIC = "fu_mic"
    BLOOD_TO_PLASMA_RATIO = "blood_to_plasma_ratio"


class AdmeStatus(StrEnum):
    MEASURED = "measured"
    NOT_MEASURED = "not_measured"
    FAILED_QC = "failed_qc"


_UNIT_TABLES: dict[AdmeAssayKind, dict[str, float]] = {
    AdmeAssayKind.PERMEABILITY: PERMEABILITY_TO_CM_S,
    AdmeAssayKind.MICROSOMAL_CLEARANCE: MICROSOMAL_TO_UL_MIN_MG,
    AdmeAssayKind.HEPATOCYTE_CLEARANCE: HEPATOCYTE_TO_UL_MIN_MCELL,
    AdmeAssayKind.FU_PLASMA: UNBOUND_FRACTION_TO_FRACTION,
    AdmeAssayKind.FU_MIC: UNBOUND_FRACTION_TO_FRACTION,
    AdmeAssayKind.BLOOD_TO_PLASMA_RATIO: BLOOD_TO_PLASMA_TO_RATIO,
}
CANONICAL_ADME_UNIT: dict[AdmeAssayKind, str] = {
    AdmeAssayKind.PERMEABILITY: "cm/s",
    AdmeAssayKind.MICROSOMAL_CLEARANCE: "uL/min/mg_protein",
    AdmeAssayKind.HEPATOCYTE_CLEARANCE: "uL/min/million_cells",
    AdmeAssayKind.FU_PLASMA: "fraction",
    AdmeAssayKind.FU_MIC: "fraction",
    AdmeAssayKind.BLOOD_TO_PLASMA_RATIO: "ratio",
}
_UNBOUND_FRACTIONS = {AdmeAssayKind.FU_PLASMA, AdmeAssayKind.FU_MIC}


class InVitroAdmeRecord(Record):
    """One in vitro ADME result. ``evidence_grade`` stays observed."""

    record_id: Identifier
    compound_id: Identifier
    assay_kind: AdmeAssayKind
    system: str = Field(min_length=1)
    status: AdmeStatus
    limitations: list[str] = Field(min_length=1)
    value: FiniteFloat | None = None
    unit: str | None = None
    replicate_id: str | None = None
    qc_reason: str | None = None
    evidence_grade: Literal["observed"] = "observed"
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def status_value_and_unit(self) -> "InVitroAdmeRecord":
        if RESEARCH_USE_LIMITATION not in self.limitations:
            raise ValueError(f"limitations must include {RESEARCH_USE_LIMITATION!r}.")
        if self.replicate_id is not None and not self.replicate_id.strip():
            raise ValueError("replicate_id cannot be blank.")
        if self.status is AdmeStatus.MEASURED:
            if self.value is None or self.unit is None:
                raise ValueError("A measured in vitro result requires value and unit.")
            if self.value <= 0:
                raise ValueError("A measured in vitro result must be positive.")
            if self.qc_reason is not None:
                raise ValueError("qc_reason is recorded only when status is failed_qc.")
            _require_vitro_unit(self.assay_kind, self.unit)
            if self.assay_kind in _UNBOUND_FRACTIONS:
                canonical = self.value * UNBOUND_FRACTION_TO_FRACTION[self.unit]
                if canonical > 1.0:
                    raise ValueError(
                        f"{self.assay_kind.value} must be an unbound fraction in (0, 1]; "
                        f"{self.value} {self.unit} is {canonical:g} as a fraction."
                    )
            return self
        if self.value is not None or self.unit is not None:
            raise ValueError("not_measured and failed_qc do not carry a numeric result.")
        if self.status is AdmeStatus.FAILED_QC and (self.qc_reason is None or not self.qc_reason.strip()):
            raise ValueError("failed_qc requires qc_reason.")
        if self.status is AdmeStatus.NOT_MEASURED and self.qc_reason is not None:
            raise ValueError("not_measured does not carry qc_reason.")
        return self

    def canonical_value(self) -> tuple[float, str] | None:
        """Return ``(value, unit)`` in the canonical unit of the assay, or ``None`` when not measured."""

        if self.status is not AdmeStatus.MEASURED or self.value is None or self.unit is None:
            return None
        factor = _UNIT_TABLES[self.assay_kind][self.unit]
        return self.value * factor, CANONICAL_ADME_UNIT[self.assay_kind]


class InVitroAdmeSet(Record):
    """Machine-readable DMPK contract exported to ``schemas/dmpk_records.schema.json``.

    ``physiology_profiles`` are caller-supplied assumptions, and
    ``ivive_predictions`` are exploratory predicted records that cite both an
    in vitro record and a profile in this set.
    """

    records: list[InVitroAdmeRecord] = Field(default_factory=list)
    physiology_profiles: list[PhysiologicalScalingParams] = Field(default_factory=list)
    ivive_predictions: list[IviveClearancePrediction] = Field(default_factory=list)
    schema_name: Literal["pharmacology-dmpk-records"] = "pharmacology-dmpk-records"
    schema_version: str = SCHEMA_VERSION

    @model_validator(mode="after")
    def unique_and_versioned(self) -> "InVitroAdmeSet":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        validate_adme_records(self.records)
        profile_ids = [item.profile_id for item in self.physiology_profiles]
        if len(profile_ids) != len(set(profile_ids)):
            raise ValueError("physiology profile_id values must be unique.")
        prediction_ids = [item.prediction_id for item in self.ivive_predictions]
        if len(prediction_ids) != len(set(prediction_ids)):
            raise ValueError("ivive prediction_id values must be unique.")
        record_ids = {item.record_id for item in self.records}
        for prediction in self.ivive_predictions:
            if prediction.source_record_id not in record_ids:
                raise ValueError(
                    f"{prediction.prediction_id} cites source_record_id {prediction.source_record_id}, "
                    "which is not in records."
                )
            if prediction.physiology_profile_id not in profile_ids:
                raise ValueError(
                    f"{prediction.prediction_id} cites physiology_profile_id {prediction.physiology_profile_id}, "
                    "which is not in physiology_profiles."
                )
        return self


def validate_adme_records(records: Sequence[InVitroAdmeRecord]) -> None:
    """Reject duplicate record ids. Constructing a record already checks its unit."""

    identifiers = [record.record_id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("record_id values must be unique.")


def scale_to_in_vivo_clearance(record: InVitroAdmeRecord) -> None:
    """Refuse to turn a single in vitro result into a plasma clearance.

    An exploratory hepatic clearance prediction needs explicit physiology,
    ``fu_plasma``, ``fu_inc``, and a blood-to-plasma ratio; use
    ``pharmacoweave.dmpk.ivive.predict_hepatic_clearance_well_stirred``. Its
    output stays ``predicted``.
    """

    raise AnalysisError(
        "In vitro ADME values are not scaled to an in vivo clearance from a single record. "
        f"{record.assay_kind.value} stays on the in vitro record. Use "
        "predict_hepatic_clearance_well_stirred with explicit physiology, fu_plasma, fu_inc, and "
        "blood-to-plasma ratio for an exploratory, predicted-grade IVIVE."
    )


def _require_vitro_unit(kind: AdmeAssayKind, unit: str) -> None:
    if _is_plasma_clearance(unit):
        raise ValueError(
            f"{unit} is a plasma clearance (volume/time). "
            "Microsomal clearance uses uL/min/mg_protein, hepatocyte clearance uses "
            "uL/min/million_cells, and neither is scaled to L/h."
        )
    allowed = _UNIT_TABLES[kind]
    if unit not in allowed:
        joined = ", ".join(sorted(allowed))
        raise ValueError(f"{kind.value} unit must be one of: {joined}.")


def _is_plasma_clearance(unit: str) -> bool:
    try:
        return parse_unit(unit).dimensions == CLEARANCE
    except UnitError:
        return False
