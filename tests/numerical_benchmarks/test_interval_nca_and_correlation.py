"""Equal-interval NCA and caller-supplied log-parameter correlation."""

import math
from datetime import UTC, datetime, timedelta

import pytest
from pydantic import ValidationError

from pharmacoweave.core.entities import (
    Analyte,
    AnalyteRole,
    Biosample,
    Compound,
    DataClass,
    DoseEvent,
    Observation,
    ObservationStatus,
    ObservationType,
    Route,
    StudyMetadata,
    StudyType,
    Timebase,
)
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint, equal_interval_profiles
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import iv_bolus_concentration
from pharmacoweave.pk_nca.calculate import nca_interval
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.simulation.exposure import ExposureSimulationRequest, _canonicalize, _perturb, simulate_exposure


def _config(dosing: str, route: str = "iv_bolus") -> NcaConfig:
    return NcaConfig.model_validate(
        {
            "analysis": "nca",
            "design_route": route,
            "dosing": dosing,
            "sampling": "rich_serial",
            "auc_method": "linear_up_log_down",
            "lambda_z_min_points": 3,
            "lambda_z_min_adjusted_r2": 0.9,
            "auc_inf_max_extrapolated_fraction": 0.2,
            "blq_rule": "keep_as_censored",
        }
    )


def _bolus_profile(times: list[float], concentrations: list[float], *, dosing: str, tau: float) -> ConcentrationProfile:
    points = tuple(
        QuantifiedPoint(observation_id=f"OBS-{index}", time_h=time_h, concentration_mg_l=concentration)
        for index, (time_h, concentration) in enumerate(zip(times, concentrations, strict=True))
    )
    return ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="ANALYTE-01",
        dose_event_id="DOSE-LAST",
        route=Route.IV_BOLUS,
        dose_mg=0.25,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
        tau_h=tau,
        n_doses=4,
        interval_regimen=dosing,  # type: ignore[arg-type]
    )


def test_steady_state_iv_bolus_recovers_clearance_from_auc_tau() -> None:
    tau = 2.0
    kel = 1.0
    c_peak = 5.0 / (1.0 - math.exp(-kel * tau))
    times = [0.0, 0.5, 1.0, 2.0]
    concentrations = [c_peak * math.exp(-kel * time_h) for time_h in times]
    result = nca_interval(
        _bolus_profile(times, concentrations, dosing="steady_state", tau=tau),
        _config("steady_state"),
    )
    assert result.auc_inf is None
    assert result.auc_tau == pytest.approx(5.0, rel=1e-12)
    assert result.clearance == pytest.approx(0.05, rel=1e-12)
    assert result.clearance_semantic is not None
    assert result.clearance_semantic.value == "CL"
    assert result.lambda_z == pytest.approx(1.0, rel=1e-9)
    assert result.volume == pytest.approx(0.05, rel=1e-9)
    assert result.cmin == pytest.approx(concentrations[-1], rel=1e-12)
    assert result.cavg == pytest.approx(5.0 / tau, rel=1e-12)
    assert result.fluctuation == pytest.approx((concentrations[0] - concentrations[-1]) / result.cavg)


def test_repeated_interval_does_not_report_clearance() -> None:
    tau = 2.0
    times = [0.0, 0.5, 1.0, 2.0]
    concentrations = [iv_bolus_concentration(
        dose=0.25,
        dose_unit="mg",
        volume=0.05,
        volume_unit="L",
        clearance=0.05,
        clearance_unit="L/h",
        time=time_h,
        time_unit="h",
        concentration_unit="mg/L",
    ) for time_h in times]
    result = nca_interval(
        _bolus_profile(times, concentrations, dosing="repeated_equal", tau=tau),
        _config("repeated_equal"),
    )
    assert result.clearance is None
    assert result.auc_tau == pytest.approx(5.0 * (1.0 - math.exp(-tau)), rel=1e-12)
    assert result.auc_inf is None
    assert any("only at steady state" in note for note in result.warnings)


def test_open_interval_omits_auc_tau() -> None:
    result = nca_interval(
        _bolus_profile(
            [0.0, 0.5, 1.0],
            [5.0, 5.0 * math.exp(-0.5), 5.0 * math.exp(-1.0)],
            dosing="steady_state",
            tau=2.0,
        ),
        _config("steady_state"),
    )
    assert result.auc_tau is None
    assert result.clearance is None
    assert result.auc_last is not None
    assert any("missing at 0 h or at tau" in note for note in result.warnings)


def test_equal_interval_builder_rejects_unequal_doses_and_spacing() -> None:
    study = _two_dose_study(second_dose_mg=0.5, second_dose_hours=2.0)
    with pytest.raises(AnalysisError, match="not equal"):
        equal_interval_profiles(study, regimen="repeated_equal")
    spaced = _two_dose_study(second_dose_mg=0.25, second_dose_hours=5.0, third_dose_hours=8.0)
    with pytest.raises(AnalysisError, match="intervals are not equal"):
        equal_interval_profiles(spaced, regimen="steady_state")


def test_correlation_is_reproducible_and_matches_the_log_scale() -> None:
    request = _iv_request(
        parameter_cvs={"clearance": 0.25, "volume": 0.25},
        parameter_correlation={"volume|clearance": 0.8},
        n_draws=1500,
        seed=11,
        threshold_concentration=None,
        threshold_unit=None,
    )
    first = simulate_exposure(request)
    second = simulate_exposure(request)
    assert first.cmax_p05 == pytest.approx(second.cmax_p05)
    assert first.auc_inf_p05 == pytest.approx(second.auc_inf_p05)
    state = _canonicalize(request)
    clearances: list[float] = []
    volumes: list[float] = []
    for index in range(request.n_draws):
        drawn = _perturb(state, request.parameter_cvs, request.seed, index, request.parameter_correlation)
        clearances.append(math.log(drawn.cl_l_h or 1.0))
        volumes.append(math.log(drawn.v_l or 1.0))
    assert _pearson(clearances, volumes) == pytest.approx(0.8, abs=0.05)
    assert "log-parameter" in " ".join(first.assumptions)
    assert first.evidence_type == "predicted"


def test_correlation_outside_range_or_not_positive_definite_is_rejected() -> None:
    with pytest.raises(ValidationError, match="\\(-1, 1\\)"):
        _iv_request(
            parameter_cvs={"clearance": 0.2, "volume": 0.2},
            parameter_correlation={"clearance|volume": 1.0},
            n_draws=40,
            seed=1,
        )
    with pytest.raises(AnalysisError, match="positive definite"):
        simulate_exposure(
            ExposureSimulationRequest.model_validate(
                {
                    "simulation_id": "SIM-CORR",
                    "context_of_use": "Research check of a caller-supplied log-parameter correlation.",
                    "model_id": "one_compartment_extravascular",
                    "dose": 1.0,
                    "dose_unit": "mg",
                    "ka": 1.0,
                    "ka_unit": "1/h",
                    "apparent_clearance": 1.0,
                    "apparent_clearance_unit": "L/h",
                    "apparent_volume": 1.0,
                    "apparent_volume_unit": "L",
                    "times_h": [0, 1, 4],
                    "parameter_cvs": {"ka": 0.2, "apparent_clearance": 0.2, "apparent_volume": 0.2},
                    "parameter_correlation": {
                        "ka|apparent_clearance": 0.99,
                        "ka|apparent_volume": 0.99,
                    },
                    "n_draws": 40,
                    "seed": 2,
                }
            )
        )


def _iv_request(**overrides: object) -> ExposureSimulationRequest:
    payload: dict[str, object] = {
        "simulation_id": "SIM-CORR-IV",
        "context_of_use": "Research check of a caller-supplied log-parameter correlation.",
        "model_id": "one_compartment_iv_bolus",
        "dose": 0.25,
        "dose_unit": "mg",
        "clearance": 0.05,
        "clearance_unit": "L/h",
        "volume": 0.05,
        "volume_unit": "L",
        "times_h": [0, 1, 2],
    }
    payload.update(overrides)
    return ExposureSimulationRequest.model_validate(payload)


def _pearson(left: list[float], right: list[float]) -> float:
    n = len(left)
    mean_left = sum(left) / n
    mean_right = sum(right) / n
    cov = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right, strict=True))
    var_left = sum((a - mean_left) ** 2 for a in left)
    var_right = sum((b - mean_right) ** 2 for b in right)
    return cov / math.sqrt(var_left * var_right)


def _two_dose_study(
    *,
    second_dose_mg: float,
    second_dose_hours: float,
    third_dose_hours: float | None = None,
) -> PharmacologyStudy:
    start = datetime(2024, 1, 15, 9, 0, tzinfo=UTC)
    metadata = StudyMetadata(
        study_id="PW-INTERVAL",
        protocol_id="SYNTH-INTERVAL",
        title="Synthetic equal-interval plasma PK",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Rattus norvegicus",
        design_version="1",
        description="Synthetic equal-interval concentrations. Not an animal experiment.",
        timebases=[Timebase(timebase_id="STUDY-CLOCK", description="Synthetic UTC clock.")],
    )
    doses = [
        DoseEvent(
            dose_event_id="DOSE-1",
            subject_id="SUBJ-01",
            compound_id="CMPD-1",
            absolute_time=start,
            dose_value=0.25,
            dose_unit="mg",
            route=Route.IV_BOLUS,
        ),
        DoseEvent(
            dose_event_id="DOSE-2",
            subject_id="SUBJ-01",
            compound_id="CMPD-1",
            absolute_time=start + timedelta(hours=second_dose_hours),
            dose_value=second_dose_mg,
            dose_unit="mg",
            route=Route.IV_BOLUS,
        ),
    ]
    if third_dose_hours is not None:
        doses.append(
            DoseEvent(
                dose_event_id="DOSE-3",
                subject_id="SUBJ-01",
                compound_id="CMPD-1",
                absolute_time=start + timedelta(hours=third_dose_hours),
                dose_value=0.25,
                dose_unit="mg",
                route=Route.IV_BOLUS,
            )
        )
    sample = Biosample(
        sample_id="SAMPLE-1",
        subject_id="SUBJ-01",
        collection_time_absolute=start,
        time_since_dose=0,
        time_since_dose_unit="h",
        dose_event_id="DOSE-1",
        matrix="plasma",
        timebase_id="STUDY-CLOCK",
    )
    observation = Observation(
        observation_id="OBS-1",
        sample_id="SAMPLE-1",
        observation_type=ObservationType.CONCENTRATION,
        status=ObservationStatus.QUANTIFIED,
        analyte_id="ANALYTE-1",
        assay_id="ASSAY-1",
        numeric_value=5.0,
        unit="mg/L",
    )
    return PharmacologyStudy(
        metadata=metadata,
        compounds=[Compound(compound_id="CMPD-1", name="synthetic")],
        analytes=[
            Analyte(
                analyte_id="ANALYTE-1",
                name="parent",
                analyte_role=AnalyteRole.PARENT,
                matrix="plasma",
                compound_id="CMPD-1",
            )
        ],
        dose_events=doses,
        biosamples=[sample],
        observations=[observation],
    )
