"""In vitro ADME records.

Permeability, microsomal clearance, and hepatocyte clearance are stored as
reported in vitro results. A plasma clearance unit is rejected. This module
does not scale any of them to an in vivo CL.
"""

from collections.abc import Sequence
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.units import CLEARANCE, parse_unit
from pharmacoweave.errors import AnalysisError, UnitError

PERMEABILITY_TO_CM_S = {"cm/s": 1.0, "um/s": 1e-4, "nm/s": 1e-7}
MICROSOMAL_TO_UL_MIN_MG = {"uL/min/mg_protein": 1.0, "mL/min/mg_protein": 1000.0}
HEPATOCYTE_TO_UL_MIN_MCELL = {"uL/min/million_cells": 1.0, "mL/min/million_cells": 1000.0}


class AdmeAssayKind(StrEnum):
    PERMEABILITY = "permeability"
    MICROSOMAL_CLEARANCE = "microsomal_clearance"
    HEPATOCYTE_CLEARANCE = "hepatocyte_clearance"


class AdmeStatus(StrEnum):
    MEASURED = "measured"
    NOT_MEASURED = "not_measured"
    FAILED_QC = "failed_qc"


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
            return self
        if self.value is not None or self.unit is not None:
            raise ValueError("not_measured and failed_qc do not carry a numeric result.")
        if self.status is AdmeStatus.FAILED_QC and (self.qc_reason is None or not self.qc_reason.strip()):
            raise ValueError("failed_qc requires qc_reason.")
        if self.status is AdmeStatus.NOT_MEASURED and self.qc_reason is not None:
            raise ValueError("not_measured does not carry qc_reason.")
        return self


def validate_adme_records(records: Sequence[InVitroAdmeRecord]) -> None:
    """Reject duplicate record ids. Constructing a record already checks its unit."""

    identifiers = [record.record_id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("record_id values must be unique.")


def scale_to_in_vivo_clearance(record: InVitroAdmeRecord) -> None:
    """Refuse to turn an in vitro result into a plasma clearance."""

    raise AnalysisError(
        "In vitro ADME values are not scaled to an in vivo clearance in this version. "
        f"{record.assay_kind.value} stays on the in vitro record."
    )


def _require_vitro_unit(kind: AdmeAssayKind, unit: str) -> None:
    if _is_plasma_clearance(unit):
        raise ValueError(
            f"{unit} is a plasma clearance (volume/time). "
            "Microsomal clearance uses uL/min/mg_protein, hepatocyte clearance uses "
            "uL/min/million_cells, and neither is scaled to L/h."
        )
    allowed = {
        AdmeAssayKind.PERMEABILITY: PERMEABILITY_TO_CM_S,
        AdmeAssayKind.MICROSOMAL_CLEARANCE: MICROSOMAL_TO_UL_MIN_MG,
        AdmeAssayKind.HEPATOCYTE_CLEARANCE: HEPATOCYTE_TO_UL_MIN_MCELL,
    }[kind]
    if unit not in allowed:
        joined = ", ".join(sorted(allowed))
        raise ValueError(f"{kind.value} unit must be one of: {joined}.")


def _is_plasma_clearance(unit: str) -> bool:
    try:
        return parse_unit(unit).dimensions == CLEARANCE
    except UnitError:
        return False
