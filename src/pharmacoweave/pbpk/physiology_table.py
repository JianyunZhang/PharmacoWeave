"""Caller-supplied tissue physiology. Validation stores no concentration."""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import PositiveFloat, Record
from pharmacoweave.core.units import CLEARANCE, DIMENSIONLESS, VOLUME, convert, parse_unit
from pharmacoweave.errors import AnalysisError

_FLOW_REL_TOL = 1e-6
_FLOW_ABS_TOL = 1e-9


class TissuePhysiology(Record):
    """One tissue. Unlisted tissues are omitted rather than filled with zero."""

    tissue_name: str = Field(min_length=1)
    volume: PositiveFloat
    volume_unit: str = Field(min_length=1)
    blood_flow: PositiveFloat
    blood_flow_unit: str = Field(min_length=1)
    partition_coefficient: PositiveFloat
    partition_unit: str = Field(min_length=1)
    source: str = Field(min_length=3)


class PhysiologyTable(Record):
    """Explicit physiology for one request. No species default is applied."""

    body_weight_kg: PositiveFloat
    body_weight_source: str = Field(min_length=3)
    cardiac_output: PositiveFloat
    cardiac_output_unit: str = Field(min_length=1)
    cardiac_output_source: str = Field(min_length=3)
    tissues: list[TissuePhysiology] = Field(min_length=1)


class PhysiologyTableCheck(Record):
    """A passed physiology check. This record has no concentration, AUC, or clearance."""

    table: PhysiologyTable
    passed: bool = True
    limitations: list[str]
    schema_version: Literal["0.1.0"] = "0.1.0"


def validate_physiology_table(table: PhysiologyTable) -> PhysiologyTableCheck:
    """Check units, sources, and flow balance. A passed check still has no predicted concentration."""

    names = [item.tissue_name for item in table.tissues]
    if len(names) != len(set(names)):
        raise AnalysisError("Each tissue name can appear once in a physiology table.")
    _positive_dimension(table.cardiac_output, table.cardiac_output_unit, CLEARANCE, "cardiac output")
    cardiac_l_h = convert(table.cardiac_output, table.cardiac_output_unit, "L/h")
    flows: list[float] = []
    for tissue in table.tissues:
        _positive_dimension(tissue.volume, tissue.volume_unit, VOLUME, f"{tissue.tissue_name} volume")
        _positive_dimension(tissue.blood_flow, tissue.blood_flow_unit, CLEARANCE, f"{tissue.tissue_name} blood flow")
        if parse_unit(tissue.partition_unit).dimensions != DIMENSIONLESS:
            raise AnalysisError(f"{tissue.tissue_name} partition coefficient must be dimensionless.")
        flows.append(convert(tissue.blood_flow, tissue.blood_flow_unit, "L/h"))
    total_flow = math.fsum(flows)
    if not math.isclose(total_flow, cardiac_l_h, rel_tol=_FLOW_REL_TOL, abs_tol=_FLOW_ABS_TOL):
        raise AnalysisError(
            f"Tissue blood flows sum to {total_flow} L/h, which does not match cardiac output {cardiac_l_h} L/h."
        )
    return PhysiologyTableCheck(
        table=table,
        limitations=[
            "Every physiological number was supplied by the caller with a source note.",
            "No species lookup table filled a missing tissue, volume, flow, or partition coefficient.",
            "A passed check does not produce a concentration, AUC, or clearance.",
        ],
    )


def _positive_dimension(value: float, unit: str, expected: object, label: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise AnalysisError(f"{label} must be positive.")
    if parse_unit(unit).dimensions != expected:
        raise AnalysisError(f"{label} unit {unit!r} has the wrong dimension.")
