"""Exploratory single-dose exposure intervals from tested analytical PK models.

The curve is a prediction under caller-declared parameters. It is not an NCA
result, not a fit to observations, and not a dosing recommendation. Parameter
bands, when requested, are percentiles of independent lognormal perturbations.
They are not confidence intervals, prediction intervals, or therapeutic windows.
"""

from __future__ import annotations

import math
import random
from typing import Literal, NoReturn, TypeVar

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, PositiveFloat, Record
from pharmacoweave.core.units import convert, convert_with_molecular_weight
from pharmacoweave.errors import AnalysisError, ContractViolation, PlannedCapabilityError
from pharmacoweave.pk_models.analytical import (
    bateman_mg_l,
    two_compartment_bolus_mg_l,
    two_compartment_macro_constants,
)

MODEL_IV_BOLUS = "one_compartment_iv_bolus"
MODEL_IV_INFUSION = "one_compartment_iv_infusion"
MODEL_EXTRAVASCULAR = "one_compartment_extravascular"
MODEL_TWO_COMPARTMENT = "two_compartment_iv_bolus"
AnalyticalModel = Literal[
    "one_compartment_iv_bolus",
    "one_compartment_iv_infusion",
    "one_compartment_extravascular",
    "two_compartment_iv_bolus",
]

_MIN_DRAWS = 40
_MAX_DRAWS = 5000
_MAX_TIMES = 500
_MAX_CV = 2.0
_KA_KEL_LIMIT = 1e-6

EXPOSURE_SIMULATION_LIMITATION = (
    "Exploratory single-dose simulation under the stated linear compartmental model. "
    "The curve was not estimated from observations. It is not a dosing recommendation "
    "and not a clinical exposure target."
)
_REGIMEN_LIMITATION = (
    "Exploratory linear superposition of equal, equally spaced doses. "
    "The curve was not estimated from observations. At steady state, AUC to infinity is not finite; "
    "auc_tau is the exposure over one dosing interval. It is not a dosing recommendation "
    "and not a clinical exposure target."
)
INTERVAL_MEANING = (
    "When present, p05, p50, and p95 are empirical percentiles of concentration, Cmax, "
    "or AUC_inf under independent lognormal parameter perturbations around the declared "
    "point values. They are not a confidence interval, not a prediction interval for a "
    "new subject, and not a therapeutic window."
)

_NUMERIC_PARAMETERS: dict[str, tuple[str, ...]] = {
    MODEL_IV_BOLUS: ("clearance", "volume"),
    MODEL_IV_INFUSION: ("clearance", "volume"),
    MODEL_EXTRAVASCULAR: ("ka", "apparent_clearance", "apparent_volume"),
    MODEL_TWO_COMPARTMENT: (
        "clearance",
        "central_volume",
        "intercompartmental_clearance",
        "peripheral_volume",
    ),
}
_REQUIRED_FIELDS: dict[str, tuple[str, ...]] = {
    MODEL_IV_BOLUS: ("clearance", "clearance_unit", "volume", "volume_unit"),
    MODEL_IV_INFUSION: (
        "clearance",
        "clearance_unit",
        "volume",
        "volume_unit",
        "infusion_duration",
        "infusion_duration_unit",
    ),
    MODEL_EXTRAVASCULAR: (
        "ka",
        "ka_unit",
        "apparent_clearance",
        "apparent_clearance_unit",
        "apparent_volume",
        "apparent_volume_unit",
    ),
    MODEL_TWO_COMPARTMENT: (
        "clearance",
        "clearance_unit",
        "central_volume",
        "central_volume_unit",
        "intercompartmental_clearance",
        "intercompartmental_clearance_unit",
        "peripheral_volume",
        "peripheral_volume_unit",
    ),
}
_ALL_OPTIONAL = (
    "clearance",
    "clearance_unit",
    "volume",
    "volume_unit",
    "infusion_duration",
    "infusion_duration_unit",
    "ka",
    "ka_unit",
    "apparent_clearance",
    "apparent_clearance_unit",
    "apparent_volume",
    "apparent_volume_unit",
    "central_volume",
    "central_volume_unit",
    "intercompartmental_clearance",
    "intercompartmental_clearance_unit",
    "peripheral_volume",
    "peripheral_volume_unit",
)


class DeclaredParameter(Record):
    """One declared input, kept under its pharmacokinetic role."""

    name: str
    value: float
    unit: str
    role: Literal[
        "dose",
        "clearance",
        "apparent_clearance",
        "volume",
        "apparent_volume",
        "central_volume",
        "peripheral_volume",
        "intercompartmental_clearance",
        "absorption_rate",
        "infusion_duration",
        "dosing_interval",
    ]


class SimulatedConcentration(Record):
    """Predicted concentration at one time. Percentiles are optional."""

    time_h: float
    concentration: float
    p05: float | None = None
    p50: float | None = None
    p95: float | None = None


class ExposureSimulationRequest(Record):
    """Analytical exposure simulation under an explicit context of use.

    The default is one dose. Repeated dosing and steady state require a dosing
    interval. They remain predicted curves, not dosing recommendations.
    """

    schema_name: Literal["pharmacoweave-exposure-simulation-request"] = "pharmacoweave-exposure-simulation-request"
    schema_version: Literal["0.2.1"] = "0.2.1"
    simulation_id: Identifier
    context_of_use: str = Field(min_length=16)
    model_id: AnalyticalModel
    dose: PositiveFloat
    dose_unit: str = Field(min_length=1)
    clearance: PositiveFloat | None = None
    clearance_unit: str | None = None
    volume: PositiveFloat | None = None
    volume_unit: str | None = None
    infusion_duration: PositiveFloat | None = None
    infusion_duration_unit: str | None = None
    ka: PositiveFloat | None = None
    ka_unit: str | None = None
    apparent_clearance: PositiveFloat | None = None
    apparent_clearance_unit: str | None = None
    apparent_volume: PositiveFloat | None = None
    apparent_volume_unit: str | None = None
    central_volume: PositiveFloat | None = None
    central_volume_unit: str | None = None
    intercompartmental_clearance: PositiveFloat | None = None
    intercompartmental_clearance_unit: str | None = None
    peripheral_volume: PositiveFloat | None = None
    peripheral_volume_unit: str | None = None
    times_h: list[FiniteFloat] | None = None
    time_start_h: FiniteFloat = 0.0
    time_end_h: PositiveFloat | None = None
    n_times: int = Field(default=0, ge=0, le=_MAX_TIMES)
    concentration_unit: str = "mg/L"
    auc_unit: str = "mg*h/L"
    threshold_concentration: PositiveFloat | None = None
    threshold_unit: str | None = None
    parameter_cvs: dict[str, FiniteFloat] = Field(default_factory=dict)
    n_draws: int = Field(default=0, ge=0, le=_MAX_DRAWS)
    seed: int | None = None
    molecular_weight: PositiveFloat | None = None
    molecular_weight_unit: str | None = None
    n_doses: int = Field(default=1, ge=1, le=500)
    dosing_interval: PositiveFloat | None = None
    dosing_interval_unit: str | None = None
    steady_state: bool = False

    @model_validator(mode="after")
    def lock_regimen_scope(self) -> ExposureSimulationRequest:
        if not self.context_of_use.strip():
            raise ValueError("context_of_use must describe the research question.")
        required = _REQUIRED_FIELDS[self.model_id]
        missing = [name for name in required if getattr(self, name) is None]
        if missing:
            raise ValueError(f"{self.model_id} requires {', '.join(missing)}.")
        unused = [name for name in _ALL_OPTIONAL if name not in required and getattr(self, name) is not None]
        if unused:
            raise ValueError(
                f"{self.model_id} does not accept {', '.join(unused)}. "
                "Absolute CL and V cannot be declared for an extravascular simulation."
            )
        allowed_cvs = set(_NUMERIC_PARAMETERS[self.model_id])
        unknown = sorted(set(self.parameter_cvs) - allowed_cvs)
        if unknown:
            raise ValueError(f"parameter_cvs keys {unknown} are not parameters of {self.model_id}.")
        for name, cv in self.parameter_cvs.items():
            if not 0.0 < cv <= _MAX_CV:
                raise ValueError(f"parameter_cvs[{name}] must be in (0, {_MAX_CV}].")
        if self.n_draws == 0 and self.parameter_cvs:
            raise ValueError("parameter_cvs requires n_draws >= 40 and a seed.")
        if 0 < self.n_draws < _MIN_DRAWS:
            raise ValueError(f"n_draws must be 0 or at least {_MIN_DRAWS}.")
        if self.n_draws > 0 and self.seed is None:
            raise ValueError("seed is required when n_draws > 0.")
        if self.n_draws > 0 and not self.parameter_cvs:
            raise ValueError("n_draws requires at least one parameter_cvs entry.")
        if (self.molecular_weight is None) != (self.molecular_weight_unit is None):
            raise ValueError("molecular_weight and molecular_weight_unit must be supplied together.")
        if (self.threshold_concentration is None) != (self.threshold_unit is None):
            raise ValueError("threshold_concentration and threshold_unit must be supplied together.")
        if self.times_h is not None:
            if self.n_times != 0 or self.time_end_h is not None:
                raise ValueError("times_h cannot be combined with n_times or time_end_h.")
            if len(self.times_h) < 2 or len(self.times_h) > _MAX_TIMES:
                raise ValueError(f"times_h must contain 2 to {_MAX_TIMES} points.")
            if self.times_h[0] < 0 or any(
                self.times_h[index] >= self.times_h[index + 1] for index in range(len(self.times_h) - 1)
            ):
                raise ValueError("times_h must be strictly increasing and start at or after 0.")
        else:
            if self.n_times < 2 or self.time_end_h is None:
                raise ValueError("Provide times_h, or n_times >= 2 with time_end_h.")
            if self.time_start_h < 0 or self.time_end_h <= self.time_start_h:
                raise ValueError("time_start_h must be >= 0 and time_end_h must be greater than time_start_h.")
        if (self.dosing_interval is None) != (self.dosing_interval_unit is None):
            raise ValueError("dosing_interval and dosing_interval_unit must be supplied together.")
        if self.steady_state and self.n_doses != 1:
            raise ValueError("steady_state is the limiting interval, so n_doses must stay 1.")
        if self.n_doses > 1 or self.steady_state:
            if self.dosing_interval is None:
                raise ValueError("Repeated dosing and steady state require dosing_interval and dosing_interval_unit.")
        elif self.dosing_interval is not None:
            raise ValueError("A single dose does not take a dosing_interval.")
        return self


class ExposureSimulationResult(Record):
    """Predicted exposure summary. ``evidence_type`` stays predicted."""

    schema_name: Literal["pharmacoweave-exposure-simulation"] = "pharmacoweave-exposure-simulation"
    schema_version: Literal["0.2.1"] = "0.2.1"
    simulation_id: Identifier
    model_id: AnalyticalModel
    context_of_use: str
    evidence_type: Literal["predicted"] = "predicted"
    constitutes_dosing_recommendation: Literal[False] = False
    concentration_unit: str
    auc_unit: str
    parameters: list[DeclaredParameter]
    curve: list[SimulatedConcentration]
    cmax: float
    cmax_p05: float | None = None
    cmax_p50: float | None = None
    cmax_p95: float | None = None
    tmax_h: float
    auc_inf: float | None = None
    auc_inf_p05: float | None = None
    auc_inf_p50: float | None = None
    auc_inf_p95: float | None = None
    auc_last: float
    regimen: Literal["single_dose", "repeated_dose", "steady_state"] = "single_dose"
    n_doses: int = 1
    dosing_interval_h: float | None = None
    auc_tau: float | None = None
    auc_tau_p05: float | None = None
    auc_tau_p50: float | None = None
    auc_tau_p95: float | None = None
    cmin: float | None = None
    accumulation_ratio: float | None = None
    half_life_h: float
    terminal_rate_name: Literal["kel", "ka", "ka_equals_kel", "beta"]
    time_above_threshold_h: float | None = None
    threshold_concentration: float | None = None
    threshold_unit: str | None = None
    n_draws: int = 0
    seed: int | None = None
    interval_meaning: str = INTERVAL_MEANING
    assumptions: list[str]
    limitations: list[str]


class _State:
    """Canonical mg / L / h parameters for one draw."""

    def __init__(
        self,
        *,
        model_id: str,
        dose_mg: float,
        cl_l_h: float | None = None,
        v_l: float | None = None,
        infusion_h: float | None = None,
        ka_h: float | None = None,
        clf_l_h: float | None = None,
        vf_l: float | None = None,
        q_l_h: float | None = None,
        v2_l: float | None = None,
    ) -> None:
        self.model_id = model_id
        self.dose_mg = dose_mg
        self.cl_l_h = cl_l_h
        self.v_l = v_l
        self.infusion_h = infusion_h
        self.ka_h = ka_h
        self.clf_l_h = clf_l_h
        self.vf_l = vf_l
        self.q_l_h = q_l_h
        self.v2_l = v2_l

    def copy(self) -> _State:
        return _State(
            model_id=self.model_id,
            dose_mg=self.dose_mg,
            cl_l_h=self.cl_l_h,
            v_l=self.v_l,
            infusion_h=self.infusion_h,
            ka_h=self.ka_h,
            clf_l_h=self.clf_l_h,
            vf_l=self.vf_l,
            q_l_h=self.q_l_h,
            v2_l=self.v2_l,
        )


class _Regimen:
    """Equal, equally spaced doses. ``tau_h`` is None only for one dose."""

    def __init__(self, kind: Literal["single", "repeated", "steady"], n_doses: int, tau_h: float | None) -> None:
        self.kind = kind
        self.n_doses = n_doses
        self.tau_h = tau_h


def simulate_exposure(request: ExposureSimulationRequest) -> ExposureSimulationResult:
    """Simulate an explicit linear concentration curve and its closed-form summaries."""

    state = _canonicalize(request)
    regimen = _regimen_from_request(request, state)
    times = _time_grid(request)
    _check_regimen_times(times, regimen)
    point_curve = [_regimen_concentration(state, time_h, regimen) for time_h in times]
    bands = _uncertainty_bands(request, state, times, regimen)
    cmax, tmax = _regimen_cmax_tmax(state, regimen)
    single_auc = _auc_inf_mg_h_l(state)
    auc_inf, auc_tau = _regimen_auc(single_auc, regimen)
    half_life, terminal_name = _terminal_half_life(state)
    threshold_mg_l = _threshold_mg_l(request)
    time_above = None if threshold_mg_l is None else _regimen_time_above(state, regimen, threshold_mg_l)
    cmin = _regimen_cmin(state, regimen)
    single_cmax, _single_tmax = _cmax_tmax(state)
    accumulation = None if regimen.kind == "single" or single_cmax <= 0 else cmax / single_cmax

    def cmax_metric(item: _State) -> float:
        return _regimen_cmax_tmax(item, regimen)[0]

    cmax_band = None if bands is None else _summary_band(request, state, cmax_metric)
    auc_source = None
    if bands is not None:
        auc_source = _summary_band(request, state, lambda item: _auc_inf_mg_h_l(item))
    concentration_unit = request.concentration_unit
    auc_unit = request.auc_unit
    weight = request.molecular_weight
    weight_unit = request.molecular_weight_unit
    curve: list[SimulatedConcentration] = []
    for index, time_h in enumerate(times):
        low = mid = high = None
        if bands is not None:
            low, mid, high = (_to_concentration(item, concentration_unit, weight, weight_unit) for item in bands[index])
        curve.append(
            SimulatedConcentration(
                time_h=time_h,
                concentration=_to_concentration(point_curve[index], concentration_unit, weight, weight_unit),
                p05=low,
                p50=mid,
                p95=high,
            )
        )
    cmax_out = _to_concentration(cmax, concentration_unit, weight, weight_unit)
    cmin_out = None if cmin is None else _to_concentration(cmin, concentration_unit, weight, weight_unit)
    auc_out = None if auc_inf is None else _to_auc(auc_inf, auc_unit, weight, weight_unit)
    auc_tau_out = None if auc_tau is None else _to_auc(auc_tau, auc_unit, weight, weight_unit)
    auc_last_out = _to_auc(_auc_last_mg_h_l(times, point_curve), auc_unit, weight, weight_unit)
    cmax_p05 = cmax_p50 = cmax_p95 = None
    auc_p05 = auc_p50 = auc_p95 = None
    if cmax_band is not None:
        cmax_p05 = _to_concentration(cmax_band[0], concentration_unit, weight, weight_unit)
        cmax_p50 = _to_concentration(cmax_band[1], concentration_unit, weight, weight_unit)
        cmax_p95 = _to_concentration(cmax_band[2], concentration_unit, weight, weight_unit)
    if auc_source is not None:
        scale = 0.0 if regimen.kind == "steady" else float(regimen.n_doses if regimen.kind == "repeated" else 1)
        if regimen.kind != "steady":
            auc_p05 = _to_auc(scale * auc_source[0], auc_unit, weight, weight_unit)
            auc_p50 = _to_auc(scale * auc_source[1], auc_unit, weight, weight_unit)
            auc_p95 = _to_auc(scale * auc_source[2], auc_unit, weight, weight_unit)
        else:
            auc_p05 = _to_auc(auc_source[0], auc_unit, weight, weight_unit)
            auc_p50 = _to_auc(auc_source[1], auc_unit, weight, weight_unit)
            auc_p95 = _to_auc(auc_source[2], auc_unit, weight, weight_unit)
    regimen_name: Literal["single_dose", "repeated_dose", "steady_state"]
    if regimen.kind == "repeated":
        regimen_name = "repeated_dose"
    elif regimen.kind == "steady":
        regimen_name = "steady_state"
    else:
        regimen_name = "single_dose"
    limitation = EXPOSURE_SIMULATION_LIMITATION if regimen.kind == "single" else _REGIMEN_LIMITATION
    return ExposureSimulationResult(
        simulation_id=request.simulation_id,
        model_id=request.model_id,
        context_of_use=request.context_of_use.strip(),
        concentration_unit=concentration_unit,
        auc_unit=auc_unit,
        parameters=_declared_parameters(request),
        curve=curve,
        cmax=cmax_out,
        cmax_p05=cmax_p05,
        cmax_p50=cmax_p50,
        cmax_p95=cmax_p95,
        tmax_h=tmax,
        auc_inf=auc_out,
        auc_inf_p05=None if regimen.kind == "steady" else auc_p05,
        auc_inf_p50=None if regimen.kind == "steady" else auc_p50,
        auc_inf_p95=None if regimen.kind == "steady" else auc_p95,
        auc_last=auc_last_out,
        regimen=regimen_name,
        n_doses=1 if regimen.kind == "steady" else regimen.n_doses,
        dosing_interval_h=regimen.tau_h,
        auc_tau=auc_tau_out,
        auc_tau_p05=auc_p05 if regimen.kind == "steady" else None,
        auc_tau_p50=auc_p50 if regimen.kind == "steady" else None,
        auc_tau_p95=auc_p95 if regimen.kind == "steady" else None,
        cmin=cmin_out,
        accumulation_ratio=accumulation,
        half_life_h=half_life,
        terminal_rate_name=terminal_name,
        time_above_threshold_h=time_above,
        threshold_concentration=request.threshold_concentration,
        threshold_unit=request.threshold_unit,
        n_draws=request.n_draws,
        seed=request.seed,
        assumptions=_assumptions(request, state, regimen),
        limitations=[limitation, INTERVAL_MEANING],
    )


def dosing_recommendation(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse to turn a simulated exposure into a dose."""

    raise ContractViolation(
        "An exploratory exposure simulation is not a dosing recommendation and cannot be converted into one."
    )


def human_starting_dose_from_simulation(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse a first-in-human or regimen dose derived from this simulation."""

    raise ContractViolation(
        "human_starting_dose_from_simulation is refused. The exposure interval is exploratory and "
        "does not justify a human starting dose."
    )


def refuse_unscoped_design(model_id: str) -> NoReturn:
    """Refuse population, PBPK, and QSP simulation."""

    raise PlannedCapabilityError(
        module="simulation",
        phase="P4",
        detail=(
            f"{model_id} is outside the analytical exposure simulator. "
            "Population covariance, PBPK, and QSP are not calculated here."
        ),
    )


def _canonicalize(request: ExposureSimulationRequest) -> _State:
    dose_mg = convert_with_molecular_weight(
        request.dose,
        request.dose_unit,
        "mg",
        molecular_weight=request.molecular_weight,
        molecular_weight_unit=request.molecular_weight_unit,
    )
    if request.model_id == MODEL_IV_BOLUS:
        return _State(
            model_id=request.model_id,
            dose_mg=dose_mg,
            cl_l_h=convert(_required(request.clearance, "clearance"), _required(request.clearance_unit, "unit"), "L/h"),
            v_l=convert(_required(request.volume, "volume"), _required(request.volume_unit, "unit"), "L"),
        )
    if request.model_id == MODEL_IV_INFUSION:
        return _State(
            model_id=request.model_id,
            dose_mg=dose_mg,
            cl_l_h=convert(_required(request.clearance, "clearance"), _required(request.clearance_unit, "unit"), "L/h"),
            v_l=convert(_required(request.volume, "volume"), _required(request.volume_unit, "unit"), "L"),
            infusion_h=convert(
                _required(request.infusion_duration, "infusion_duration"),
                _required(request.infusion_duration_unit, "unit"),
                "h",
            ),
        )
    if request.model_id == MODEL_EXTRAVASCULAR:
        return _State(
            model_id=request.model_id,
            dose_mg=dose_mg,
            ka_h=convert(_required(request.ka, "ka"), _required(request.ka_unit, "unit"), "1/h"),
            clf_l_h=convert(
                _required(request.apparent_clearance, "apparent_clearance"),
                _required(request.apparent_clearance_unit, "unit"),
                "L/h",
            ),
            vf_l=convert(
                _required(request.apparent_volume, "apparent_volume"),
                _required(request.apparent_volume_unit, "unit"),
                "L",
            ),
        )
    return _State(
        model_id=request.model_id,
        dose_mg=dose_mg,
        cl_l_h=convert(_required(request.clearance, "clearance"), _required(request.clearance_unit, "unit"), "L/h"),
        v_l=convert(
            _required(request.central_volume, "central_volume"),
            _required(request.central_volume_unit, "unit"),
            "L",
        ),
        q_l_h=convert(
            _required(request.intercompartmental_clearance, "intercompartmental_clearance"),
            _required(request.intercompartmental_clearance_unit, "unit"),
            "L/h",
        ),
        v2_l=convert(
            _required(request.peripheral_volume, "peripheral_volume"),
            _required(request.peripheral_volume_unit, "unit"),
            "L",
        ),
    )


def _declared_parameters(request: ExposureSimulationRequest) -> list[DeclaredParameter]:
    rows = [
        DeclaredParameter(name="dose", value=request.dose, unit=request.dose_unit, role="dose"),
    ]
    if request.model_id in {MODEL_IV_BOLUS, MODEL_IV_INFUSION}:
        rows.append(
            DeclaredParameter(
                name="CL",
                value=_required(request.clearance, "clearance"),
                unit=_required(request.clearance_unit, "unit"),
                role="clearance",
            )
        )
        rows.append(
            DeclaredParameter(
                name="V",
                value=_required(request.volume, "volume"),
                unit=_required(request.volume_unit, "unit"),
                role="volume",
            )
        )
    if request.model_id == MODEL_IV_INFUSION:
        rows.append(
            DeclaredParameter(
                name="infusion_duration",
                value=_required(request.infusion_duration, "infusion_duration"),
                unit=_required(request.infusion_duration_unit, "unit"),
                role="infusion_duration",
            )
        )
    if request.model_id == MODEL_EXTRAVASCULAR:
        rows.extend(
            [
                DeclaredParameter(
                    name="ka",
                    value=_required(request.ka, "ka"),
                    unit=_required(request.ka_unit, "unit"),
                    role="absorption_rate",
                ),
                DeclaredParameter(
                    name="CL_over_F",
                    value=_required(request.apparent_clearance, "apparent_clearance"),
                    unit=_required(request.apparent_clearance_unit, "unit"),
                    role="apparent_clearance",
                ),
                DeclaredParameter(
                    name="V_over_F",
                    value=_required(request.apparent_volume, "apparent_volume"),
                    unit=_required(request.apparent_volume_unit, "unit"),
                    role="apparent_volume",
                ),
            ]
        )
    if request.model_id == MODEL_TWO_COMPARTMENT:
        rows.extend(
            [
                DeclaredParameter(
                    name="CL",
                    value=_required(request.clearance, "clearance"),
                    unit=_required(request.clearance_unit, "unit"),
                    role="clearance",
                ),
                DeclaredParameter(
                    name="V1",
                    value=_required(request.central_volume, "central_volume"),
                    unit=_required(request.central_volume_unit, "unit"),
                    role="central_volume",
                ),
                DeclaredParameter(
                    name="Q",
                    value=_required(request.intercompartmental_clearance, "Q"),
                    unit=_required(request.intercompartmental_clearance_unit, "unit"),
                    role="intercompartmental_clearance",
                ),
                DeclaredParameter(
                    name="V2",
                    value=_required(request.peripheral_volume, "peripheral_volume"),
                    unit=_required(request.peripheral_volume_unit, "unit"),
                    role="peripheral_volume",
                ),
            ]
        )
    if request.dosing_interval is not None and request.dosing_interval_unit is not None:
        rows.append(
            DeclaredParameter(
                name="tau",
                value=request.dosing_interval,
                unit=request.dosing_interval_unit,
                role="dosing_interval",
            )
        )
    return rows


def _assumptions(request: ExposureSimulationRequest, state: _State, regimen: _Regimen) -> list[str]:
    if regimen.kind == "single":
        regimen_note = "Linear kinetics and a single dose."
    elif regimen.kind == "steady":
        regimen_note = (
            f"Linear kinetics at steady state. Doses are equal and spaced by {regimen.tau_h:g} h. "
            "Times are measured from the start of one dosing interval."
        )
    else:
        regimen_note = (
            f"Linear superposition of {regimen.n_doses} equal doses spaced by {regimen.tau_h:g} h. "
            "No inter-dose covariance was applied."
        )
    notes = [
        regimen_note,
        "The curve uses the declared parameters. It was not fitted to study observations.",
        "Bioavailability, protein binding, and interindividual covariance were not inferred.",
    ]
    if request.model_id == MODEL_EXTRAVASCULAR:
        notes.append("Bioavailability is unknown, so clearance and volume are CL/F and V/F, not CL and V.")
        if state.ka_h is not None and state.clf_l_h is not None and state.vf_l is not None:
            kel = state.clf_l_h / state.vf_l
            if abs(state.ka_h - kel) <= _KA_KEL_LIMIT * max(state.ka_h, kel):
                notes.append("ka equals kel within the numerical limit, so the Bateman limit form was used.")
            elif state.ka_h < kel:
                notes.append("ka is slower than kel, so the terminal half-life is governed by ka (flip-flop).")
    if request.model_id == MODEL_IV_INFUSION:
        notes.append(
            "Input is a zero-order infusion of the declared duration. The infused amount is the declared dose."
        )
    if request.model_id == MODEL_TWO_COMPARTMENT:
        notes.append("The reported concentration is the central compartment. Peripheral concentration is not reported.")
    if request.n_draws == 0:
        notes.append("No parameter uncertainty was propagated.")
    else:
        varied = ", ".join(f"{name} CV={cv:g}" for name, cv in sorted(request.parameter_cvs.items()))
        notes.append(
            f"{request.n_draws} independent lognormal draws (seed {request.seed}) perturbed {varied}. "
            "Parameters without a CV stayed fixed. Draws are independent, so this is not a covariance matrix."
        )
    return notes


def _time_grid(request: ExposureSimulationRequest) -> list[float]:
    if request.times_h is not None:
        return list(request.times_h)
    assert request.time_end_h is not None
    step = (request.time_end_h - request.time_start_h) / (request.n_times - 1)
    return [request.time_start_h + index * step for index in range(request.n_times)]


def _concentration_mg_l(state: _State, time_h: float) -> float:
    if time_h < 0:
        raise AnalysisError("Simulation time must be >= 0.")
    if state.model_id == MODEL_IV_BOLUS:
        kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
        return (state.dose_mg / _required(state.v_l, "V")) * math.exp(-kel * time_h)
    if state.model_id == MODEL_IV_INFUSION:
        return _infusion_mg_l(state, time_h)
    if state.model_id == MODEL_EXTRAVASCULAR:
        return bateman_mg_l(
            dose_mg=state.dose_mg,
            ka_h=_required(state.ka_h, "ka"),
            clf_l_h=_required(state.clf_l_h, "CL/F"),
            vf_l=_required(state.vf_l, "V/F"),
            time_h=time_h,
        )
    return two_compartment_bolus_mg_l(
        dose_mg=state.dose_mg,
        cl_l_h=_required(state.cl_l_h, "CL"),
        v1_l=_required(state.v_l, "V1"),
        q_l_h=_required(state.q_l_h, "Q"),
        v2_l=_required(state.v2_l, "V2"),
        time_h=time_h,
    )


def _regimen_from_request(request: ExposureSimulationRequest, state: _State) -> _Regimen:
    if request.n_doses == 1 and not request.steady_state:
        return _Regimen("single", 1, None)
    tau_h = convert(
        _required(request.dosing_interval, "dosing_interval"),
        _required(request.dosing_interval_unit, "dosing_interval_unit"),
        "h",
    )
    if state.infusion_h is not None and state.infusion_h >= tau_h:
        raise AnalysisError(
            "Infusion duration must be shorter than the dosing interval. Continuous infusion is a different design."
        )
    if request.steady_state:
        return _Regimen("steady", 1, tau_h)
    return _Regimen("repeated", request.n_doses, tau_h)


def _check_regimen_times(times: list[float], regimen: _Regimen) -> None:
    if regimen.kind != "steady" or regimen.tau_h is None:
        return
    if any(time_h < 0.0 or time_h > regimen.tau_h + 1e-9 for time_h in times):
        raise AnalysisError(
            "Steady-state times are measured from the start of one dosing interval and must lie in [0, tau]."
        )


def _regimen_concentration(state: _State, time_h: float, regimen: _Regimen) -> float:
    if regimen.kind == "single" or regimen.tau_h is None:
        return _concentration_mg_l(state, time_h)
    if regimen.kind == "steady":
        return _steady_concentration(state, min(time_h, regimen.tau_h), regimen.tau_h)
    total = 0.0
    last_index = min(regimen.n_doses - 1, int(math.floor(time_h / regimen.tau_h + 1e-10)))
    for index in range(last_index + 1):
        elapsed = time_h - index * regimen.tau_h
        if elapsed >= 0.0:
            total += _concentration_mg_l(state, elapsed)
    return total


def _regimen_auc(single_auc: float, regimen: _Regimen) -> tuple[float | None, float | None]:
    if regimen.kind == "steady":
        return None, single_auc
    if regimen.kind == "repeated":
        return regimen.n_doses * single_auc, None
    return single_auc, None


def _regimen_cmin(state: _State, regimen: _Regimen) -> float | None:
    if regimen.kind == "single" or regimen.tau_h is None:
        return None
    if regimen.kind == "steady":
        if state.model_id == MODEL_IV_INFUSION:
            return _infusion_ss(state, 0.0, regimen.tau_h)
        return _steady_concentration(state, regimen.tau_h, regimen.tau_h)
    return _regimen_concentration(state, regimen.n_doses * regimen.tau_h, regimen)


def _regimen_cmax_tmax(state: _State, regimen: _Regimen) -> tuple[float, float]:
    if regimen.kind == "single" or regimen.tau_h is None:
        return _cmax_tmax(state)
    tau = regimen.tau_h
    if state.model_id == MODEL_IV_BOLUS:
        kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
        amount_over_volume = state.dose_mg / _required(state.v_l, "V")
        if regimen.kind == "steady":
            return amount_over_volume / (1.0 - math.exp(-kel * tau)), 0.0
        cmax = amount_over_volume * (1.0 - math.exp(-kel * regimen.n_doses * tau)) / (1.0 - math.exp(-kel * tau))
        return cmax, (regimen.n_doses - 1) * tau
    if state.model_id == MODEL_IV_INFUSION:
        duration = _required(state.infusion_h, "infusion_duration")
        if regimen.kind == "steady":
            return _infusion_ss(state, duration, tau), duration
        peak_time = (regimen.n_doses - 1) * tau + duration
        return _regimen_concentration(state, peak_time, regimen), peak_time
    if state.model_id == MODEL_TWO_COMPARTMENT:
        if regimen.kind == "steady":
            return _two_comp_ss(state, 0.0, tau), 0.0
        peak_time = (regimen.n_doses - 1) * tau
        return _regimen_concentration(state, peak_time, regimen), peak_time
    if regimen.kind == "steady":
        tmax = _oral_ss_tmax(state, tau)
        return _oral_ss(state, tmax, tau), tmax
    left = (regimen.n_doses - 1) * tau
    return _maximize(lambda time_h: _regimen_concentration(state, time_h, regimen), left, regimen.n_doses * tau)


def _regimen_time_above(state: _State, regimen: _Regimen, threshold_mg_l: float) -> float:
    if regimen.kind == "single" or regimen.tau_h is None:
        return _time_above_h(state, threshold_mg_l)
    if state.model_id == MODEL_IV_BOLUS:
        return _bolus_regimen_time_above(state, regimen, threshold_mg_l)
    horizon = regimen.tau_h if regimen.kind == "steady" else regimen.n_doses * regimen.tau_h
    return _time_above_scan(state, regimen, threshold_mg_l, horizon)


def _bolus_regimen_time_above(state: _State, regimen: _Regimen, threshold_mg_l: float) -> float:
    tau = _required(regimen.tau_h, "dosing_interval")
    kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
    amount_over_volume = state.dose_mg / _required(state.v_l, "V")
    if regimen.kind == "steady":
        cmax = amount_over_volume / (1.0 - math.exp(-kel * tau))
        cmin = cmax * math.exp(-kel * tau)
        if threshold_mg_l >= cmax:
            return 0.0
        if threshold_mg_l <= cmin:
            return tau
        return math.log(cmax / threshold_mg_l) / kel
    total = 0.0
    for index in range(regimen.n_doses):
        cpeak = amount_over_volume * (1.0 - math.exp(-kel * (index + 1) * tau)) / (1.0 - math.exp(-kel * tau))
        if cpeak > threshold_mg_l:
            total += min(tau, math.log(cpeak / threshold_mg_l) / kel)
    return total


def _time_above_scan(state: _State, regimen: _Regimen, threshold_mg_l: float, horizon: float) -> float:
    steps = 4000
    total = 0.0
    previous_time = 0.0
    previous = _regimen_concentration(state, 0.0, regimen)
    for index in range(1, steps + 1):
        time_h = horizon * index / steps
        current = _regimen_concentration(state, time_h, regimen)
        span = time_h - previous_time
        if previous >= threshold_mg_l and current >= threshold_mg_l:
            total += span
        elif previous >= threshold_mg_l and current < threshold_mg_l and previous != current:
            total += (previous - threshold_mg_l) / (previous - current) * span
        elif previous < threshold_mg_l and current >= threshold_mg_l and previous != current:
            total += (current - threshold_mg_l) / (current - previous) * span
        previous_time = time_h
        previous = current
    return total


def _steady_concentration(state: _State, time_h: float, tau_h: float) -> float:
    if state.model_id == MODEL_IV_BOLUS:
        kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
        amount_over_volume = state.dose_mg / _required(state.v_l, "V")
        return amount_over_volume * math.exp(-kel * time_h) / (1.0 - math.exp(-kel * tau_h))
    if state.model_id == MODEL_IV_INFUSION:
        return _infusion_ss(state, time_h, tau_h)
    if state.model_id == MODEL_EXTRAVASCULAR:
        return _oral_ss(state, time_h, tau_h)
    return _two_comp_ss(state, time_h, tau_h)


def _infusion_ss(state: _State, time_h: float, tau_h: float) -> float:
    duration = _required(state.infusion_h, "infusion_duration")
    clearance = _required(state.cl_l_h, "CL")
    volume = _required(state.v_l, "V")
    kel = clearance / volume
    rate = state.dose_mg / duration
    cmax = (rate / clearance) * (1.0 - math.exp(-kel * duration)) / (1.0 - math.exp(-kel * tau_h))
    cmin = cmax * math.exp(-kel * (tau_h - duration))
    if time_h <= duration:
        return (rate / clearance) * -math.expm1(-kel * time_h) + cmin * math.exp(-kel * time_h)
    return cmax * math.exp(-kel * (time_h - duration))


def _oral_ss(state: _State, time_h: float, tau_h: float) -> float:
    ka = _required(state.ka_h, "ka")
    kel = _required(state.clf_l_h, "CL/F") / _required(state.vf_l, "V/F")
    dose_over_volume = state.dose_mg / _required(state.vf_l, "V/F")
    if abs(ka - kel) <= _KA_KEL_LIMIT * max(ka, kel):
        rate = ka
        ratio = math.exp(-rate * tau_h)
        geometric = time_h / (1.0 - ratio) + tau_h * ratio / (1.0 - ratio) ** 2
        return dose_over_volume * rate * math.exp(-rate * time_h) * geometric
    factor = ka / (ka - kel) * dose_over_volume
    return factor * (
        math.exp(-kel * time_h) / (1.0 - math.exp(-kel * tau_h))
        - math.exp(-ka * time_h) / (1.0 - math.exp(-ka * tau_h))
    )


def _oral_ss_tmax(state: _State, tau_h: float) -> float:
    ka = _required(state.ka_h, "ka")
    kel = _required(state.clf_l_h, "CL/F") / _required(state.vf_l, "V/F")
    if abs(ka - kel) <= _KA_KEL_LIMIT * max(ka, kel):
        _peak, tmax = _maximize(lambda time_h: _oral_ss(state, time_h, tau_h), 0.0, tau_h)
        return tmax
    numerator = ka / (1.0 - math.exp(-ka * tau_h))
    denominator = kel / (1.0 - math.exp(-kel * tau_h))
    tmax = math.log(numerator / denominator) / (ka - kel)
    if 0.0 < tmax < tau_h:
        return tmax
    _peak, tmax = _maximize(lambda time_h: _oral_ss(state, time_h, tau_h), 0.0, tau_h)
    return tmax


def _two_comp_ss(state: _State, time_h: float, tau_h: float) -> float:
    coefficient_a, alpha, coefficient_b, beta = two_compartment_macro_constants(
        dose_mg=state.dose_mg,
        cl_l_h=_required(state.cl_l_h, "CL"),
        v1_l=_required(state.v_l, "V1"),
        q_l_h=_required(state.q_l_h, "Q"),
        v2_l=_required(state.v2_l, "V2"),
    )
    return coefficient_a * math.exp(-alpha * time_h) / (1.0 - math.exp(-alpha * tau_h)) + coefficient_b * math.exp(
        -beta * time_h
    ) / (1.0 - math.exp(-beta * tau_h))


def _maximize(function, left: float, right: float) -> tuple[float, float]:
    best_time = left
    best_value = function(left)
    steps = 400
    for index in range(1, steps + 1):
        time_h = left + (right - left) * index / steps
        value = function(time_h)
        if value > best_value:
            best_value = value
            best_time = time_h
    return best_value, best_time


def _infusion_mg_l(state: _State, time_h: float) -> float:
    duration = _required(state.infusion_h, "infusion_duration")
    clearance = _required(state.cl_l_h, "CL")
    volume = _required(state.v_l, "V")
    if state.dose_mg == 0.0 or time_h == 0.0:
        return 0.0
    kel = clearance / volume
    rate = state.dose_mg / duration
    active = min(time_h, duration)
    at_active = -(rate / clearance) * math.expm1(-kel * active)
    if time_h <= duration:
        return at_active
    return at_active * math.exp(-kel * (time_h - duration))


def _cmax_tmax(state: _State) -> tuple[float, float]:
    if state.model_id == MODEL_IV_BOLUS:
        return state.dose_mg / _required(state.v_l, "V"), 0.0
    if state.model_id == MODEL_IV_INFUSION:
        duration = _required(state.infusion_h, "infusion_duration")
        return _infusion_mg_l(state, duration), duration
    if state.model_id == MODEL_TWO_COMPARTMENT:
        return state.dose_mg / _required(state.v_l, "V1"), 0.0
    ka = _required(state.ka_h, "ka")
    kel = _required(state.clf_l_h, "CL/F") / _required(state.vf_l, "V/F")
    if abs(ka - kel) <= _KA_KEL_LIMIT * max(ka, kel):
        tmax = 1.0 / ka
    else:
        tmax = math.log(ka / kel) / (ka - kel)
    return _concentration_mg_l(state, tmax), tmax


def _auc_inf_mg_h_l(state: _State) -> float:
    if state.model_id == MODEL_EXTRAVASCULAR:
        return state.dose_mg / _required(state.clf_l_h, "CL/F")
    return state.dose_mg / _required(state.cl_l_h, "CL")


def _terminal_half_life(state: _State) -> tuple[float, Literal["kel", "ka", "ka_equals_kel", "beta"]]:
    if state.model_id == MODEL_TWO_COMPARTMENT:
        _a, _alpha, _b, beta = two_compartment_macro_constants(
            dose_mg=state.dose_mg,
            cl_l_h=_required(state.cl_l_h, "CL"),
            v1_l=_required(state.v_l, "V1"),
            q_l_h=_required(state.q_l_h, "Q"),
            v2_l=_required(state.v2_l, "V2"),
        )
        return math.log(2.0) / beta, "beta"
    if state.model_id == MODEL_EXTRAVASCULAR:
        ka = _required(state.ka_h, "ka")
        kel = _required(state.clf_l_h, "CL/F") / _required(state.vf_l, "V/F")
        if abs(ka - kel) <= _KA_KEL_LIMIT * max(ka, kel):
            return math.log(2.0) / ka, "ka_equals_kel"
        if ka < kel:
            return math.log(2.0) / ka, "ka"
        return math.log(2.0) / kel, "kel"
    kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
    return math.log(2.0) / kel, "kel"


def _time_above_h(state: _State, threshold_mg_l: float) -> float:
    cmax, tmax = _cmax_tmax(state)
    if threshold_mg_l >= cmax:
        return 0.0
    if state.model_id == MODEL_IV_BOLUS:
        kel = _required(state.cl_l_h, "CL") / _required(state.v_l, "V")
        return math.log(cmax / threshold_mg_l) / kel
    if state.model_id == MODEL_IV_INFUSION:
        duration = _required(state.infusion_h, "infusion_duration")
        clearance = _required(state.cl_l_h, "CL")
        volume = _required(state.v_l, "V")
        kel = clearance / volume
        rate = state.dose_mg / duration
        remaining = 1.0 - threshold_mg_l * clearance / rate
        if remaining <= 0.0:
            return 0.0
        start = -math.log(remaining) / kel
        end = duration + math.log(cmax / threshold_mg_l) / kel
        return max(0.0, end - start)
    if state.model_id == MODEL_TWO_COMPARTMENT:
        _a, _alpha, _b, beta = two_compartment_macro_constants(
            dose_mg=state.dose_mg,
            cl_l_h=_required(state.cl_l_h, "CL"),
            v1_l=_required(state.v_l, "V1"),
            q_l_h=_required(state.q_l_h, "Q"),
            v2_l=_required(state.v2_l, "V2"),
        )
        upper = _past_threshold(state, math.log(cmax / threshold_mg_l) / beta, threshold_mg_l)
        return _bisection_descending(state, 0.0, upper, threshold_mg_l)
    start = _bisection_rising(state, 0.0, tmax, threshold_mg_l)
    kel = _required(state.clf_l_h, "CL/F") / _required(state.vf_l, "V/F")
    ka = _required(state.ka_h, "ka")
    slow = min(ka, kel)
    upper = _past_threshold(state, tmax + math.log(cmax / threshold_mg_l) / slow + 1.0, threshold_mg_l)
    end = _bisection_descending(state, tmax, upper, threshold_mg_l)
    return max(0.0, end - start)


def _past_threshold(state: _State, upper: float, threshold_mg_l: float) -> float:
    for _ in range(20):
        if _concentration_mg_l(state, upper) <= threshold_mg_l:
            return upper
        upper *= 2.0
    raise AnalysisError("The concentration did not fall through the declared threshold.")


def _bisection_rising(state: _State, left: float, right: float, target: float) -> float:
    for _ in range(80):
        mid = 0.5 * (left + right)
        if _concentration_mg_l(state, mid) < target:
            left = mid
        else:
            right = mid
    return 0.5 * (left + right)


def _bisection_descending(state: _State, left: float, right: float, target: float) -> float:
    for _ in range(80):
        mid = 0.5 * (left + right)
        if _concentration_mg_l(state, mid) > target:
            left = mid
        else:
            right = mid
    return 0.5 * (left + right)


def _threshold_mg_l(request: ExposureSimulationRequest) -> float | None:
    if request.threshold_concentration is None or request.threshold_unit is None:
        return None
    return convert_with_molecular_weight(
        request.threshold_concentration,
        request.threshold_unit,
        "mg/L",
        molecular_weight=request.molecular_weight,
        molecular_weight_unit=request.molecular_weight_unit,
    )


def _uncertainty_bands(
    request: ExposureSimulationRequest,
    state: _State,
    times: list[float],
    regimen: _Regimen,
) -> list[tuple[float, float, float]] | None:
    if request.n_draws == 0:
        return None
    draws = [_perturb(state, request.parameter_cvs, request.seed, index) for index in range(request.n_draws)]
    bands: list[tuple[float, float, float]] = []
    for time_index in range(len(times)):
        samples = [_regimen_concentration(draw, times[time_index], regimen) for draw in draws]
        bands.append(_three_percentiles(samples))
    return bands


def _summary_band(request: ExposureSimulationRequest, state: _State, metric) -> tuple[float, float, float]:
    draws = [_perturb(state, request.parameter_cvs, request.seed, index) for index in range(request.n_draws)]
    return _three_percentiles([metric(draw) for draw in draws])


def _perturb(state: _State, cvs: dict[str, float], seed: int | None, index: int) -> _State:
    generator = random.Random((0 if seed is None else seed) * 1_000_003 + index + 1)
    drawn = state.copy()
    mapping = {
        "clearance": "cl_l_h",
        "volume": "v_l",
        "central_volume": "v_l",
        "peripheral_volume": "v2_l",
        "intercompartmental_clearance": "q_l_h",
        "ka": "ka_h",
        "apparent_clearance": "clf_l_h",
        "apparent_volume": "vf_l",
    }
    for name, cv in cvs.items():
        attribute = mapping[name]
        point = getattr(drawn, attribute)
        sigma = math.sqrt(math.log1p(cv * cv))
        setattr(drawn, attribute, point * math.exp(sigma * generator.gauss(0.0, 1.0)))
    return drawn


def _three_percentiles(values: list[float]) -> tuple[float, float, float]:
    return _percentile(values, 0.05), _percentile(values, 0.50), _percentile(values, 0.95)


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    weight = position - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def _auc_last_mg_h_l(times: list[float], concentrations: list[float]) -> float:
    total = 0.0
    for index in range(len(times) - 1):
        delta = times[index + 1] - times[index]
        total += 0.5 * delta * (concentrations[index] + concentrations[index + 1])
    return total


def _to_concentration(
    value_mg_l: float,
    unit: str,
    molecular_weight: float | None,
    molecular_weight_unit: str | None,
) -> float:
    return convert_with_molecular_weight(
        value_mg_l,
        "mg/L",
        unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def _to_auc(
    value_mg_h_l: float,
    unit: str,
    molecular_weight: float | None,
    molecular_weight_unit: str | None,
) -> float:
    return convert_with_molecular_weight(
        value_mg_h_l,
        "mg*h/L",
        unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


_T = TypeVar("_T")


def _required(value: _T | None, name: str) -> _T:
    if value is None:
        raise AnalysisError(f"{name} is required.")
    return value
