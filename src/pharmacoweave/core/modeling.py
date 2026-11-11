"""Model specifications and result records.

These records describe an analysis. Constructing one does not fit a model.
A passing diagnostic must carry limitations, so goodness of fit alone cannot
mark a result as acceptable. IV clearance and extravascular CL/F stay on
different semantic labels.
"""

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Identifier, Record, Route
from pharmacoweave.core.units import (
    AUC_AMOUNT,
    AUC_MASS,
    CLEARANCE,
    CONC_AMOUNT,
    CONC_MASS,
    RATE,
    TIME,
    VOLUME,
    parse_unit,
)
from pharmacoweave.version import SCHEMA_VERSION


class ParameterSemantic(StrEnum):
    CL = "CL"
    CL_OVER_F = "CL_over_F"
    V = "V"
    V_OVER_F = "V_over_F"
    VZ = "Vz"
    VZ_OVER_F = "Vz_over_F"
    V1 = "V1"
    V2 = "V2"
    V3 = "V3"
    Q = "Q"
    Q2 = "Q2"
    Q3 = "Q3"
    KA = "ka"
    KEL = "kel"
    K12 = "k12"
    K21 = "k21"
    HALF_LIFE = "half_life"
    F = "F"


class ErrorModel(StrEnum):
    ADDITIVE = "additive"
    PROPORTIONAL = "proportional"
    COMBINED = "combined"
    UNSPECIFIED = "unspecified"


class DiagnosticStatus(StrEnum):
    NOT_RUN = "not_run"
    PASS = "pass"
    WARNING = "warning"
    FAIL = "fail"


_IV_ROUTES = {Route.IV_BOLUS, Route.IV_INFUSION}
_EXTRAVASCULAR_ROUTES = {Route.ORAL, Route.SUBCUTANEOUS, Route.INTRAMUSCULAR}
_ABSOLUTE_CLEARANCE = {
    ParameterSemantic.CL,
    ParameterSemantic.V,
    ParameterSemantic.VZ,
    ParameterSemantic.V1,
    ParameterSemantic.V2,
    ParameterSemantic.V3,
    ParameterSemantic.Q,
    ParameterSemantic.Q2,
    ParameterSemantic.Q3,
}
_APPARENT_CLEARANCE = {
    ParameterSemantic.CL_OVER_F,
    ParameterSemantic.V_OVER_F,
    ParameterSemantic.VZ_OVER_F,
}
_SEMANTIC_DIMENSIONS = {
    ParameterSemantic.CL: CLEARANCE,
    ParameterSemantic.CL_OVER_F: CLEARANCE,
    ParameterSemantic.V: VOLUME,
    ParameterSemantic.V_OVER_F: VOLUME,
    ParameterSemantic.VZ: VOLUME,
    ParameterSemantic.VZ_OVER_F: VOLUME,
    ParameterSemantic.V1: VOLUME,
    ParameterSemantic.V2: VOLUME,
    ParameterSemantic.V3: VOLUME,
    ParameterSemantic.Q: CLEARANCE,
    ParameterSemantic.Q2: CLEARANCE,
    ParameterSemantic.Q3: CLEARANCE,
    ParameterSemantic.KA: RATE,
    ParameterSemantic.KEL: RATE,
    ParameterSemantic.K12: RATE,
    ParameterSemantic.K21: RATE,
    ParameterSemantic.HALF_LIFE: TIME,
    ParameterSemantic.F: (0, 0, 0, 0, 0),
}


class ModelSpec(Record):
    """A named structural model and the question it is allowed to answer."""

    model_id: Identifier
    model_type: str = Field(min_length=1)
    equations_ref: str = Field(min_length=1)
    parameters: list[str] = Field(min_length=1)
    error_model: ErrorModel
    assumptions: list[str] = Field(min_length=1)
    question_of_interest: str = Field(min_length=10)
    context_of_use: str = Field(min_length=10)
    estimation_implemented: bool = False


class Metric(Record):
    """One named diagnostic number."""

    name: str = Field(min_length=1)
    value: float
    unit: str | None = None


class DiagnosticRecord(Record):
    """A model diagnostic. Pass requires limitations to be written down."""

    fit_id: Identifier
    status: DiagnosticStatus
    limitations: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    metrics: list[Metric] = Field(default_factory=list)

    @model_validator(mode="after")
    def pass_states_limits(self) -> "DiagnosticRecord":
        if self.status == DiagnosticStatus.PASS and not self.limitations:
            raise ValueError("A passing diagnostic must state limitations. Goodness of fit alone is not acceptance.")
        return self


class ModelFit(Record):
    """One fit against one dataset. Fits are appended, not overwritten."""

    fit_id: Identifier
    model_id: Identifier
    dataset_id: Identifier
    software_version: str = Field(min_length=1)
    estimation_method: str = Field(min_length=1)
    context_of_use: str = Field(min_length=10)
    diagnostic_status: DiagnosticStatus = DiagnosticStatus.NOT_RUN
    converged: bool | None = None
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def pass_requires_convergence_and_limits(self) -> "ModelFit":
        if self.diagnostic_status == DiagnosticStatus.PASS:
            if self.converged is not True:
                raise ValueError("A fit marked pass must have converged.")
            if not self.limitations:
                raise ValueError("A fit marked pass must record limitations.")
        return self


class ParameterEstimate(Record):
    """A parameter from a fit, including the route that fixes its meaning."""

    fit_id: Identifier
    parameter: str = Field(min_length=1)
    semantic: ParameterSemantic
    value: float
    unit: str
    estimation_method: str = Field(min_length=1)
    route: Route
    standard_error: float | None = Field(default=None, ge=0)
    ci_lower: float | None = None
    ci_upper: float | None = None
    ci_level: float | None = Field(default=None, gt=0, lt=1)
    bound_hit: bool = False

    @model_validator(mode="after")
    def route_and_unit(self) -> "ParameterEstimate":
        _reject_clearance_label(self.route, self.semantic)
        expected = _SEMANTIC_DIMENSIONS[self.semantic]
        if parse_unit(self.unit).dimensions != expected:
            raise ValueError(f"Unit {self.unit} does not match {self.semantic.value}.")
        if self.semantic == ParameterSemantic.F and not 0 < self.value <= 1:
            raise ValueError("F must be in (0, 1].")
        if self.ci_lower is not None and self.ci_upper is not None and self.ci_lower > self.ci_upper:
            raise ValueError("ci_lower must be less than or equal to ci_upper.")
        if (self.ci_lower is None) != (self.ci_upper is None):
            raise ValueError("ci_lower and ci_upper are recorded together.")
        if (self.ci_lower is not None) and self.ci_level is None:
            raise ValueError("A confidence interval requires ci_level.")
        return self


class NcaParameterSet(Record):
    """NCA parameters that are allowed to be reported for one subject profile."""

    subject_id: Identifier
    analyte_id: Identifier
    route: Route
    dose_event_id: Identifier
    auc_method: str = Field(min_length=1)
    terminal_phase_qualified: bool
    auc_last: float | None = None
    auc_last_unit: str | None = None
    auc_inf: float | None = None
    auc_inf_unit: str | None = None
    auc_inf_extrap_fraction: float | None = Field(default=None, ge=0, le=1)
    cmax: float | None = None
    cmax_unit: str | None = None
    tmax: float | None = None
    tmax_unit: str | None = None
    lambda_z: float | None = None
    lambda_z_unit: str | None = None
    half_life: float | None = None
    half_life_unit: str | None = None
    terminal_phase_n_points: int | None = Field(default=None, ge=0)
    terminal_phase_adjusted_r2: float | None = Field(default=None, ge=0, le=1)
    clearance_semantic: ParameterSemantic | None = None
    clearance: float | None = None
    clearance_unit: str | None = None
    volume_semantic: ParameterSemantic | None = None
    volume: float | None = None
    volume_unit: str | None = None
    regimen: Literal["single_dose", "repeated_equal", "steady_state"] = "single_dose"
    n_doses: int | None = Field(default=None, ge=1)
    dosing_interval_h: float | None = Field(default=None, gt=0)
    auc_tau: float | None = None
    auc_tau_unit: str | None = None
    cmin: float | None = None
    cmin_unit: str | None = None
    cavg: float | None = None
    cavg_unit: str | None = None
    fluctuation: float | None = None
    accumulation_ratio: float | None = Field(default=None, gt=0)
    excluded_reason: str | None = None
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def reporting_rules(self) -> "NcaParameterSet":
        _pair(self.auc_last, self.auc_last_unit, "auc_last", {AUC_MASS, AUC_AMOUNT})
        _pair(self.cmax, self.cmax_unit, "cmax", {CONC_MASS, CONC_AMOUNT})
        if (self.tmax is None) != (self.tmax_unit is None):
            raise ValueError("tmax requires both a value and a unit.")
        if self.tmax is not None and self.tmax < 0:
            raise ValueError("tmax must be >= 0 when it is reported.")
        if self.tmax_unit is not None and parse_unit(self.tmax_unit).dimensions not in {TIME}:
            raise ValueError(f"tmax unit {self.tmax_unit} has the wrong dimension.")
        if not self.terminal_phase_qualified:
            leaked = [
                name
                for name, value in (
                    ("auc_inf", self.auc_inf),
                    ("lambda_z", self.lambda_z),
                    ("half_life", self.half_life),
                )
                if value is not None
            ]
            if leaked:
                joined = ", ".join(leaked)
                raise ValueError(f"Do not report {joined} when the terminal phase is not qualified.")
        _pair(self.auc_inf, self.auc_inf_unit, "auc_inf", {AUC_MASS, AUC_AMOUNT})
        _pair(self.lambda_z, self.lambda_z_unit, "lambda_z", {RATE})
        _pair(self.half_life, self.half_life_unit, "half_life", {TIME})
        _pair(self.auc_tau, self.auc_tau_unit, "auc_tau", {AUC_MASS, AUC_AMOUNT})
        _pair(self.cmin, self.cmin_unit, "cmin", {CONC_MASS, CONC_AMOUNT})
        _pair(self.cavg, self.cavg_unit, "cavg", {CONC_MASS, CONC_AMOUNT})
        if self.regimen == "single_dose" and any(
            value is not None
            for value in (self.auc_tau, self.cmin, self.cavg, self.fluctuation, self.dosing_interval_h, self.n_doses)
        ):
            raise ValueError("Interval NCA fields are only recorded for repeated_equal or steady_state.")
        if self.regimen == "repeated_equal" and self.clearance is not None:
            raise ValueError("Dose/AUC_tau clearance is recorded only when the interval is declared steady_state.")
        if self.auc_inf is not None and self.auc_inf_extrap_fraction is None:
            raise ValueError("auc_inf requires auc_inf_extrap_fraction.")
        if self.clearance is not None or self.clearance_semantic is not None or self.clearance_unit is not None:
            if self.clearance is None or self.clearance_semantic is None or self.clearance_unit is None:
                raise ValueError("clearance, clearance_semantic, and clearance_unit are recorded together.")
            _reject_clearance_label(self.route, self.clearance_semantic)
            if self.clearance_semantic not in {ParameterSemantic.CL, ParameterSemantic.CL_OVER_F}:
                raise ValueError("NCA clearance_semantic must be CL or CL_over_F.")
            if parse_unit(self.clearance_unit).dimensions != CLEARANCE:
                raise ValueError("clearance_unit must be a clearance.")
        if self.volume is not None or self.volume_semantic is not None or self.volume_unit is not None:
            if self.volume is None or self.volume_semantic is None or self.volume_unit is None:
                raise ValueError("volume, volume_semantic, and volume_unit are recorded together.")
            _reject_clearance_label(self.route, self.volume_semantic)
            allowed_volume = {
                ParameterSemantic.V,
                ParameterSemantic.V_OVER_F,
                ParameterSemantic.VZ,
                ParameterSemantic.VZ_OVER_F,
            }
            if self.volume_semantic not in allowed_volume:
                raise ValueError("NCA volume_semantic must be V, V_over_F, Vz, or Vz_over_F.")
            if parse_unit(self.volume_unit).dimensions != VOLUME:
                raise ValueError("volume_unit must be a volume.")
        return self


class AnalysisRecords(Record):
    """Machine-readable analysis records. An empty document means no analysis ran."""

    model_specs: list[ModelSpec] = Field(default_factory=list)
    model_fits: list[ModelFit] = Field(default_factory=list)
    parameter_estimates: list[ParameterEstimate] = Field(default_factory=list)
    nca_parameters: list[NcaParameterSet] = Field(default_factory=list)
    diagnostics: list[DiagnosticRecord] = Field(default_factory=list)
    schema_name: str = "pharmacology-analysis-records"
    schema_version: str = SCHEMA_VERSION

    @model_validator(mode="after")
    def contract_version(self) -> "AnalysisRecords":
        if self.schema_name != "pharmacology-analysis-records":
            raise ValueError("schema_name must be pharmacology-analysis-records.")
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        return self


def auc_inf_is_reportable(record: NcaParameterSet, *, max_extrapolated_fraction: float = 0.20) -> bool:
    """Return whether AUC∞ meets the terminal-phase and extrapolation policy.

    The default 0.20 fraction is a reporting flag, not a claim that a larger
    extrapolation is biologically impossible. Callers that report AUC∞ use
    this helper instead of treating a present number as sufficient.
    """

    if not 0 < max_extrapolated_fraction < 1:
        raise ValueError("max_extrapolated_fraction must be between 0 and 1.")
    if record.auc_inf is None or record.auc_inf_extrap_fraction is None:
        return False
    return record.terminal_phase_qualified and record.auc_inf_extrap_fraction <= max_extrapolated_fraction


def _pair(value: float | None, unit: str | None, label: str, allowed: set[tuple[int, int, int, int, int]]) -> None:
    if value is None and unit is None:
        return
    if value is None or unit is None:
        raise ValueError(f"{label} requires both a value and a unit.")
    if value <= 0:
        raise ValueError(f"{label} must be positive when it is reported.")
    if parse_unit(unit).dimensions not in allowed:
        raise ValueError(f"{label} unit {unit} has the wrong dimension.")


def _reject_clearance_label(route: Route, semantic: ParameterSemantic) -> None:
    if route in _IV_ROUTES and semantic in _APPARENT_CLEARANCE:
        raise ValueError(f"Route {route.value} cannot be reported as {semantic.value}. Use CL and V.")
    if route in _EXTRAVASCULAR_ROUTES and semantic in _ABSOLUTE_CLEARANCE:
        raise ValueError(
            f"Route {route.value} cannot be reported as {semantic.value}. Use CL/F and V/F when F is unknown."
        )
    if route == Route.OTHER and semantic in _ABSOLUTE_CLEARANCE | _APPARENT_CLEARANCE:
        raise ValueError(
            f"Route other cannot carry {semantic.value}. Use a specific route before labeling clearance or volume."
        )
