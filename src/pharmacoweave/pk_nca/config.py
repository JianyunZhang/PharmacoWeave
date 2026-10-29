"""NCA configuration for the one design this version is prepared to accept."""

from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Record
from pharmacoweave.version import SCHEMA_VERSION


class NcaConfig(Record):
    """Reporting rules for single-dose IV bolus rich sampling.

    Loading a config does not calculate NCA. Other designs fail validation
    here so they are not silently passed to a later calculator.
    """

    analysis: Literal["nca"]
    design_route: Literal["iv_bolus"]
    dosing: Literal["single_dose"]
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
