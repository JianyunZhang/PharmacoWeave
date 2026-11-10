"""Steady-state algebra for liver plus one non-eliminating tissue.

The liver is well stirred with fu = 1:

    CL = Q_liver * CLint / (Q_liver + CLint)
    C_arterial = input_rate / CL
    C_liver = Kp_liver * C_arterial * Q_liver / (Q_liver + CLint)
    C_tissue = Kp_tissue * C_arterial

No differential equation is integrated.
"""

from typing import Literal

from pharmacoweave.core.entities import Record
from pharmacoweave.core.units import CLEARANCE, convert, parse_unit
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pbpk.physiology_table import PhysiologyTable, validate_physiology_table

ClearanceName = Literal["CL", "CL/F"]


class SteadyStatePbpk(Record):
    """Predicted steady-state concentrations. Physiological numbers remain assumptions."""

    evidence_grade: Literal["predicted"] = "predicted"
    clearance_name: ClearanceName
    hepatic_clearance_l_h: float
    arterial_concentration_mg_l: float
    liver_concentration_mg_l: float
    tissue_name: str
    tissue_concentration_mg_l: float
    assumption_sources: list[str]
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def steady_state_one_tissue(
    table: PhysiologyTable,
    *,
    input_rate: float,
    input_rate_unit: str,
    hepatic_intrinsic_clearance: float,
    clearance_unit: str,
    clearance_name: ClearanceName,
    time_course: bool = False,
) -> SteadyStatePbpk:
    """Calculate the algebraic steady state. A requested time course is rejected."""

    if time_course:
        raise AnalysisError("A time course is rejected. This function reports steady state only.")
    if input_rate <= 0 or hepatic_intrinsic_clearance <= 0:
        raise AnalysisError("Input rate and hepatic intrinsic clearance must be positive.")
    check = validate_physiology_table(table)
    if len(table.tissues) != 2:
        raise AnalysisError("Steady state accepts the liver and exactly one other tissue.")
    by_name = {item.tissue_name: item for item in table.tissues}
    liver = by_name.get("liver")
    if liver is None:
        raise AnalysisError("The physiology table must name one tissue 'liver'.")
    other = next(item for item in table.tissues if item.tissue_name != "liver")
    if parse_unit(input_rate_unit).dimensions != parse_unit("mg/h").dimensions:
        raise AnalysisError("input_rate_unit must be a drug-mass per time unit such as mg/h.")
    if parse_unit(clearance_unit).dimensions != CLEARANCE:
        raise AnalysisError("clearance_unit must be a clearance unit such as L/h.")
    rate_mg_h = convert(input_rate, input_rate_unit, "mg/h")
    clint = convert(hepatic_intrinsic_clearance, clearance_unit, "L/h")
    liver_flow = convert(liver.blood_flow, liver.blood_flow_unit, "L/h")
    hepatic_clearance = liver_flow * clint / (liver_flow + clint)
    arterial = rate_mg_h / hepatic_clearance
    liver_concentration = liver.partition_coefficient * arterial * liver_flow / (liver_flow + clint)
    tissue_concentration = other.partition_coefficient * arterial
    sources = [
        table.body_weight_source,
        table.cardiac_output_source,
        liver.source,
        other.source,
    ]
    return SteadyStatePbpk(
        clearance_name=clearance_name,
        hepatic_clearance_l_h=hepatic_clearance,
        arterial_concentration_mg_l=arterial,
        liver_concentration_mg_l=liver_concentration,
        tissue_name=other.tissue_name,
        tissue_concentration_mg_l=tissue_concentration,
        assumption_sources=sources,
        limitations=[
            *check.limitations,
            "fu is 1. Liver clearance is Q_liver * CLint / (Q_liver + CLint).",
            "Arterial concentration is input rate divided by that hepatic clearance.",
            "Liver concentration is Kp_liver times arterial concentration times Q_liver / (Q_liver + CLint).",
            "The other tissue concentration is Kp times arterial concentration.",
            "The physiological numbers are caller-supplied assumptions.",
            f"The clearance parameter name remains {clearance_name}.",
            "The result is predicted and is not a dosing recommendation.",
        ],
    )
