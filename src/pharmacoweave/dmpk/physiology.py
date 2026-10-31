"""Explicit physiological assumptions and exploratory IVIVE prediction records.

Nothing here contains a species lookup table. Every scaling factor that an
IVIVE uses (body weight, liver weight, hepatic blood flow, MPPGL or
hepatocellularity) has to be passed in by the caller together with a source
note. A prediction built from these values is locked to
``evidence_grade="predicted"`` and is never an observed or estimated CL.
"""

import math
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Identifier, PositiveFloat, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION

IVIVE_EXPLORATORY_LIMITATION = (
    "Exploratory IVIVE prediction from in vitro data and caller-supplied physiology. "
    "Not an observed or estimated in vivo clearance, not a starting-dose or FIH recommendation."
)

_REL_TOL = 1e-9


class HepaticClearanceModel(StrEnum):
    WELL_STIRRED = "well_stirred"


class PhysiologicalScalingParams(Record):
    """Caller-supplied physiology for one species/individual. No defaults are filled in."""

    profile_id: Identifier
    species: str = Field(min_length=1)
    body_weight_kg: PositiveFloat
    liver_weight_g: PositiveFloat
    hepatic_blood_flow_l_h: PositiveFloat
    source: str = Field(min_length=3, description="Citation or provenance of the assumed values.")
    microsomal_protein_mg_per_g_liver: PositiveFloat | None = None
    hepatocellularity_million_cells_per_g_liver: PositiveFloat | None = None
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def has_a_scaling_factor(self) -> "PhysiologicalScalingParams":
        if self.microsomal_protein_mg_per_g_liver is None and self.hepatocellularity_million_cells_per_g_liver is None:
            raise ValueError(
                "Provide microsomal_protein_mg_per_g_liver (MPPGL) or "
                "hepatocellularity_million_cells_per_g_liver (HPGL). No default is assumed."
            )
        return self


class IviveClearancePrediction(Record):
    """One exploratory hepatic clearance prediction. ``evidence_grade`` stays predicted."""

    prediction_id: Identifier
    compound_id: Identifier
    species: str = Field(min_length=1)
    physiology_profile_id: Identifier
    source_record_id: Identifier
    assay_kind: Literal["microsomal_clearance", "hepatocyte_clearance"]
    clint_app_value: PositiveFloat
    clint_app_unit: Literal["uL/min/mg_protein", "uL/min/million_cells"]
    fu_inc: PositiveFloat
    fu_plasma: PositiveFloat
    blood_to_plasma_ratio: PositiveFloat
    fu_blood: PositiveFloat
    clint_u_vitro_value: PositiveFloat
    scaling_factor_value: PositiveFloat
    scaling_factor_unit: Literal["mg_protein/g_liver", "million_cells/g_liver"]
    liver_weight_g: PositiveFloat
    body_weight_kg: PositiveFloat
    hepatic_blood_flow_l_h: PositiveFloat
    clint_u_hepatic_l_h: PositiveFloat
    hepatic_extraction_ratio: float = Field(gt=0, lt=1)
    hepatic_blood_clearance_l_h: PositiveFloat
    hepatic_plasma_clearance_l_h: PositiveFloat
    hepatic_plasma_clearance_l_h_kg: PositiveFloat
    assumptions: list[str] = Field(min_length=1)
    limitations: list[str] = Field(min_length=2)
    derived_from: list[str] = Field(min_length=1)
    model: HepaticClearanceModel = HepaticClearanceModel.WELL_STIRRED
    clearance_unit: Literal["L/h"] = "L/h"
    evidence_grade: Literal["predicted"] = "predicted"
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def equations_are_consistent(self) -> "IviveClearancePrediction":
        for name in ("fu_inc", "fu_plasma"):
            if getattr(self, name) > 1:
                raise ValueError(f"{name} must be in (0, 1].")
        for required in (RESEARCH_USE_LIMITATION, IVIVE_EXPLORATORY_LIMITATION):
            if required not in self.limitations:
                raise ValueError(f"limitations must include {required!r}.")
        if any(not item.strip() for item in self.derived_from):
            raise ValueError("derived_from entries cannot be blank.")
        checks = {
            "fu_blood": (self.fu_blood, self.fu_plasma / self.blood_to_plasma_ratio),
            "clint_u_vitro_value": (self.clint_u_vitro_value, self.clint_app_value / self.fu_inc),
            "hepatic_extraction_ratio": (
                self.hepatic_extraction_ratio,
                self.hepatic_blood_clearance_l_h / self.hepatic_blood_flow_l_h,
            ),
            "hepatic_plasma_clearance_l_h": (
                self.hepatic_plasma_clearance_l_h,
                self.blood_to_plasma_ratio * self.hepatic_blood_clearance_l_h,
            ),
            "hepatic_plasma_clearance_l_h_kg": (
                self.hepatic_plasma_clearance_l_h_kg,
                self.hepatic_plasma_clearance_l_h / self.body_weight_kg,
            ),
        }
        for name, (stored, expected) in checks.items():
            if not math.isclose(stored, expected, rel_tol=_REL_TOL):
                raise ValueError(f"{name}={stored} is inconsistent with the well-stirred equations ({expected}).")
        return self
