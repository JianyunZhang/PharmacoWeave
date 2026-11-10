"""NCA configuration for rich-serial plasma designs.

Supported routes are ``iv_bolus``, ``iv_infusion``, ``oral``,
``subcutaneous``, and ``intramuscular``. ``single_dose`` is one dose.
``repeated_equal`` and ``steady_state`` require equal doses and equal
intervals. Sparse destructive sampling, urine excretion, and ``other`` routes
fail validation here instead of being passed to a calculator.
"""

from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Record
from pharmacoweave.version import SCHEMA_VERSION


class NcaConfig(Record):
    """Reporting rules for rich-serial plasma NCA.

    ``repeated_equal`` reports the last closed interdose interval and does not
    treat Dose/AUC_tau as clearance. ``steady_state`` is a caller declaration
    that the analyzed interval is at steady state, so clearance is Dose/AUC_tau.
    """

    analysis: Literal["nca"]
    design_route: Literal["iv_bolus", "iv_infusion", "oral", "subcutaneous", "intramuscular"]
    dosing: Literal["single_dose", "repeated_equal", "steady_state"]
    sampling: Literal["rich_serial"]
    auc_method: Literal["linear_up_log_down"]
    lambda_z_min_points: int = Field(ge=3, le=20)
    lambda_z_min_adjusted_r2: float = Field(gt=0, le=1)
    auc_inf_max_extrapolated_fraction: float = Field(gt=0, lt=1)
    blq_rule: Literal["keep_as_censored", "exclude"]
    schema_version: Literal["0.1.0"] = "0.1.0"
    note: str = ""

    @model_validator(mode="after")
    def version_matches_package(self) -> "NcaConfig":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        return self
