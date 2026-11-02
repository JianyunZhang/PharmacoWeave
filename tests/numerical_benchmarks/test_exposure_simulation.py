"""Closed-form checks for exploratory single-dose exposure simulation."""

import math

import pytest
from pydantic import ValidationError

from pharmacoweave.errors import AnalysisError, ContractViolation, PlannedCapabilityError
from pharmacoweave.pk_models.analytical import bateman_mg_l, iv_infusion_concentration, two_compartment_macro_constants
from pharmacoweave.simulation.exposure import (
    ExposureSimulationRequest,
    dosing_recommendation,
    human_starting_dose_from_simulation,
    refuse_unscoped_design,
    simulate_exposure,
)


def _iv_request(**overrides: object) -> ExposureSimulationRequest:
    payload: dict[str, object] = {
        "simulation_id": "SIM-IV",
        "context_of_use": "Research check of a declared one-compartment IV bolus curve.",
        "model_id": "one_compartment_iv_bolus",
        "dose": 0.25,
        "dose_unit": "mg",
        "clearance": 0.05,
        "clearance_unit": "L/h",
        "volume": 0.05,
        "volume_unit": "L",
        "times_h": [0, 1, 2, 4],
        "threshold_concentration": 1.0,
        "threshold_unit": "mg/L",
    }
    payload.update(overrides)
    return ExposureSimulationRequest.model_validate(payload)


def test_iv_bolus_matches_dose_over_parameters() -> None:
    result = simulate_exposure(_iv_request())
    assert result.evidence_type == "predicted"
    assert result.constitutes_dosing_recommendation is False
    assert result.cmax == pytest.approx(5.0)
    assert result.tmax_h == 0.0
    assert result.auc_inf == pytest.approx(5.0)
    assert result.half_life_h == pytest.approx(math.log(2.0))
    assert result.terminal_rate_name == "kel"
    assert result.time_above_threshold_h == pytest.approx(math.log(5.0))
    assert result.curve[0].concentration == pytest.approx(5.0)
    assert result.curve[1].concentration == pytest.approx(5.0 * math.exp(-1.0))
    assert {item.name for item in result.parameters} >= {"CL", "V", "dose"}


def test_infusion_cmax_matches_the_analytical_reference() -> None:
    request = ExposureSimulationRequest.model_validate(
        {
            "simulation_id": "SIM-INF",
            "context_of_use": "Research check of a declared constant-rate IV infusion.",
            "model_id": "one_compartment_iv_infusion",
            "dose": 1.0,
            "dose_unit": "mg",
            "clearance": 0.2,
            "clearance_unit": "L/h",
            "volume": 1.0,
            "volume_unit": "L",
            "infusion_duration": 0.5,
            "infusion_duration_unit": "h",
            "times_h": [0, 0.5, 2],
        }
    )
    result = simulate_exposure(request)
    expected = iv_infusion_concentration(
        dose=1.0,
        dose_unit="mg",
        infusion_duration=0.5,
        infusion_duration_unit="h",
        volume=1.0,
        volume_unit="L",
        clearance=0.2,
        clearance_unit="L/h",
        time=0.5,
        time_unit="h",
        concentration_unit="mg/L",
    )
    assert result.cmax == pytest.approx(expected, rel=1e-12)
    assert result.tmax_h == pytest.approx(0.5)
    assert result.auc_inf == pytest.approx(5.0)
    assert result.curve[1].concentration == pytest.approx(expected, rel=1e-12)


def test_extravascular_reports_apparent_parameters_only() -> None:
    ka = 1.2
    kel = 0.2
    tmax = math.log(ka / kel) / (ka - kel)
    expected_cmax = bateman_mg_l(dose_mg=10.0, ka_h=ka, clf_l_h=2.0, vf_l=10.0, time_h=tmax)
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-EV",
                "context_of_use": "Research check of a declared extravascular curve with unknown F.",
                "model_id": "one_compartment_extravascular",
                "dose": 10.0,
                "dose_unit": "mg",
                "ka": ka,
                "ka_unit": "1/h",
                "apparent_clearance": 2.0,
                "apparent_clearance_unit": "L/h",
                "apparent_volume": 10.0,
                "apparent_volume_unit": "L",
                "n_times": 5,
                "time_end_h": 8.0,
            }
        )
    )
    names = {item.name for item in result.parameters}
    assert names == {"dose", "ka", "CL_over_F", "V_over_F"}
    assert "CL" not in names
    assert "V" not in names
    assert result.cmax == pytest.approx(expected_cmax, rel=1e-12)
    assert result.tmax_h == pytest.approx(tmax, rel=1e-12)
    assert result.auc_inf == pytest.approx(5.0)
    assert result.terminal_rate_name == "kel"


def test_flip_flop_labels_the_terminal_rate_as_ka() -> None:
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-FLIP",
                "context_of_use": "Research check that a slower ka is the terminal rate.",
                "model_id": "one_compartment_extravascular",
                "dose": 1.0,
                "dose_unit": "mg",
                "ka": 0.1,
                "ka_unit": "1/h",
                "apparent_clearance": 1.0,
                "apparent_clearance_unit": "L/h",
                "apparent_volume": 1.0,
                "apparent_volume_unit": "L",
                "times_h": [0, 1, 10],
            }
        )
    )
    assert result.terminal_rate_name == "ka"
    assert result.half_life_h == pytest.approx(math.log(2.0) / 0.1)


def test_equal_ka_and_kel_uses_the_limit_form() -> None:
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-LIMIT",
                "context_of_use": "Research check of the Bateman limit when ka equals kel.",
                "model_id": "one_compartment_extravascular",
                "dose": 1.0,
                "dose_unit": "mg",
                "ka": 0.2,
                "ka_unit": "1/h",
                "apparent_clearance": 0.2,
                "apparent_clearance_unit": "L/h",
                "apparent_volume": 1.0,
                "apparent_volume_unit": "L",
                "times_h": [0, 5, 10],
            }
        )
    )
    assert result.terminal_rate_name == "ka_equals_kel"
    assert result.tmax_h == pytest.approx(5.0)
    assert any("Bateman limit" in note for note in result.assumptions)


def test_two_compartment_cmax_and_auc_use_central_volume_and_clearance() -> None:
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-2C",
                "context_of_use": "Research check of a declared two-compartment IV bolus curve.",
                "model_id": "two_compartment_iv_bolus",
                "dose": 1.0,
                "dose_unit": "mg",
                "clearance": 0.25,
                "clearance_unit": "L/h",
                "central_volume": 0.5,
                "central_volume_unit": "L",
                "intercompartmental_clearance": 0.5,
                "intercompartmental_clearance_unit": "L/h",
                "peripheral_volume": 1.0,
                "peripheral_volume_unit": "L",
                "times_h": [0, 1, 4],
                "threshold_concentration": 0.5,
                "threshold_unit": "mg/L",
            }
        )
    )
    assert result.cmax == pytest.approx(2.0)
    assert result.auc_inf == pytest.approx(4.0)
    assert result.terminal_rate_name == "beta"
    assert result.time_above_threshold_h is not None
    assert result.time_above_threshold_h > 0.0


def test_volume_uncertainty_moves_cmax_and_not_auc() -> None:
    request = _iv_request(
        parameter_cvs={"volume": 0.2},
        n_draws=80,
        seed=7,
        threshold_concentration=None,
        threshold_unit=None,
    )
    first = simulate_exposure(request)
    second = simulate_exposure(request)
    assert first.cmax_p05 == pytest.approx(second.cmax_p05)
    assert first.cmax_p05 is not None and first.cmax_p50 is not None and first.cmax_p95 is not None
    assert first.cmax_p05 < first.cmax_p50 < first.cmax_p95
    assert first.auc_inf_p05 == pytest.approx(first.auc_inf)
    assert first.auc_inf_p50 == pytest.approx(first.auc_inf)
    assert first.auc_inf_p95 == pytest.approx(first.auc_inf)
    assert first.curve[0].p05 == pytest.approx(first.cmax_p05)
    assert "independent" in first.interval_meaning


def test_clearance_uncertainty_moves_auc_and_not_cmax() -> None:
    result = simulate_exposure(_iv_request(parameter_cvs={"clearance": 0.25}, n_draws=80, seed=11))
    assert result.cmax_p05 == pytest.approx(result.cmax)
    assert result.cmax_p95 == pytest.approx(result.cmax)
    assert result.auc_inf_p05 is not None and result.auc_inf_p95 is not None
    assert result.auc_inf_p05 < result.auc_inf < result.auc_inf_p95


def test_scope_refusals() -> None:
    with pytest.raises(ValidationError):
        _iv_request(context_of_use="too short")
    with pytest.raises(ValidationError):
        _iv_request(n_doses=2)
    with pytest.raises(ValidationError):
        _iv_request(steady_state=True)
    with pytest.raises(ValidationError, match="does not accept clearance"):
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-MIX",
                "context_of_use": "Research request that illegally mixes CL with CL/F.",
                "model_id": "one_compartment_extravascular",
                "dose": 1.0,
                "dose_unit": "mg",
                "clearance": 1.0,
                "clearance_unit": "L/h",
                "ka": 1.0,
                "ka_unit": "1/h",
                "apparent_clearance": 1.0,
                "apparent_clearance_unit": "L/h",
                "apparent_volume": 1.0,
                "apparent_volume_unit": "L",
                "times_h": [0, 1],
            }
        )
    with pytest.raises(ContractViolation):
        dosing_recommendation()
    with pytest.raises(ContractViolation):
        human_starting_dose_from_simulation()
    with pytest.raises(PlannedCapabilityError):
        refuse_unscoped_design("pbpk")


def test_bundle_persists_the_predicted_curve(tmp_path, repo_root) -> None:
    import json

    from pharmacoweave.core.study import PharmacologyStudy

    study = PharmacologyStudy.from_files(repo_root / "examples" / "data" / "iv_bolus_single_dose")
    result = simulate_exposure(_iv_request())
    destination = tmp_path / "sim-bundle"
    manifest = study.export_bundle(destination, run_id="sim-run", exposure_simulation=result)
    assert manifest.analysis_executed is True
    assert manifest.data["exposure_simulation"] == "analysis/exposure_simulation.json"
    assert not any(item.role == "exposure_simulation" for item in manifest.omitted)
    payload = json.loads((destination / "analysis" / "exposure_simulation.json").read_text(encoding="utf-8"))
    assert payload["evidence_type"] == "predicted"
    assert payload["constitutes_dosing_recommendation"] is False
    report = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "Exploratory single-dose exposure simulation" in report
    assert "not a dosing recommendation" in report


def test_repeated_iv_bolus_superposition_and_total_auc() -> None:
    result = simulate_exposure(
        _iv_request(
            n_doses=2,
            dosing_interval=1.0,
            dosing_interval_unit="h",
            times_h=[0, 1, 2],
            threshold_concentration=3.0,
            threshold_unit="mg/L",
        )
    )
    kel = 1.0
    dose_over_volume = 5.0
    cmax = dose_over_volume * (1.0 + math.exp(-kel))
    assert result.regimen == "repeated_dose"
    assert result.tmax_h == pytest.approx(1.0)
    assert result.cmax == pytest.approx(cmax)
    assert result.curve[1].concentration == pytest.approx(cmax)
    assert result.auc_inf == pytest.approx(10.0)
    assert result.auc_tau is None
    assert result.accumulation_ratio == pytest.approx(1.0 + math.exp(-1.0))
    cpeak_second = dose_over_volume * (1.0 - math.exp(-2.0 * kel)) / (1.0 - math.exp(-kel))
    expected_time = math.log(dose_over_volume / 3.0) + math.log(cpeak_second / 3.0)
    assert result.time_above_threshold_h == pytest.approx(expected_time)
    assert {item.name for item in result.parameters} >= {"CL", "V", "tau"}


def test_steady_state_iv_bolus_interval_exposure() -> None:
    tau = 2.0
    kel = 1.0
    result = simulate_exposure(
        _iv_request(
            steady_state=True,
            dosing_interval=tau,
            dosing_interval_unit="h",
            times_h=[0, tau],
            threshold_concentration=None,
            threshold_unit=None,
        )
    )
    cmax = 5.0 / (1.0 - math.exp(-kel * tau))
    assert result.regimen == "steady_state"
    assert result.auc_inf is None
    assert result.auc_tau == pytest.approx(5.0)
    assert result.cmax == pytest.approx(cmax)
    assert result.cmin == pytest.approx(cmax * math.exp(-kel * tau))
    assert result.accumulation_ratio == pytest.approx(1.0 / (1.0 - math.exp(-kel * tau)))
    assert result.curve[1].concentration == pytest.approx(result.cmin)
    assert result.constitutes_dosing_recommendation is False


def test_steady_state_infusion_matches_the_closed_form() -> None:
    tau = 8.0
    duration = 0.5
    kel = 0.2
    rate = 2.0
    clearance = 0.2
    cmax = (rate / clearance) * (1.0 - math.exp(-kel * duration)) / (1.0 - math.exp(-kel * tau))
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-INF-SS",
                "context_of_use": "Research check of a steady-state constant-rate infusion interval.",
                "model_id": "one_compartment_iv_infusion",
                "dose": 1.0,
                "dose_unit": "mg",
                "clearance": clearance,
                "clearance_unit": "L/h",
                "volume": 1.0,
                "volume_unit": "L",
                "infusion_duration": duration,
                "infusion_duration_unit": "h",
                "dosing_interval": tau,
                "dosing_interval_unit": "h",
                "steady_state": True,
                "times_h": [0, duration, tau],
            }
        )
    )
    assert result.cmax == pytest.approx(cmax, rel=1e-12)
    assert result.tmax_h == pytest.approx(duration)
    assert result.cmin == pytest.approx(cmax * math.exp(-kel * (tau - duration)), rel=1e-12)
    assert result.auc_tau == pytest.approx(5.0)
    assert result.auc_inf is None
    assert result.curve[0].concentration == pytest.approx(result.cmin)


def test_steady_state_oral_auc_and_superposition() -> None:
    tau = 12.0
    ka = 1.2
    kel = 0.2
    dose = 10.0
    volume = 10.0
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-EV-SS",
                "context_of_use": "Research check of steady-state extravascular exposure with unknown F.",
                "model_id": "one_compartment_extravascular",
                "dose": dose,
                "dose_unit": "mg",
                "ka": ka,
                "ka_unit": "1/h",
                "apparent_clearance": 2.0,
                "apparent_clearance_unit": "L/h",
                "apparent_volume": volume,
                "apparent_volume_unit": "L",
                "dosing_interval": tau,
                "dosing_interval_unit": "h",
                "steady_state": True,
                "times_h": [0, 1, tau],
            }
        )
    )
    stacked = sum(
        bateman_mg_l(dose_mg=dose, ka_h=ka, clf_l_h=2.0, vf_l=volume, time_h=1.0 + index * tau) for index in range(40)
    )
    names = {item.name for item in result.parameters}
    assert "CL" not in names
    assert "CL_over_F" in names
    assert result.auc_tau == pytest.approx(5.0)
    assert result.auc_inf is None
    assert result.curve[1].concentration == pytest.approx(stacked, rel=1e-8)
    assert result.curve[0].concentration == pytest.approx(result.curve[2].concentration, rel=1e-12)
    numerator = ka / (1.0 - math.exp(-ka * tau))
    denominator = kel / (1.0 - math.exp(-kel * tau))
    assert result.tmax_h == pytest.approx(math.log(numerator / denominator) / (ka - kel), rel=1e-12)


def test_steady_state_two_compartment_auc_and_peak() -> None:
    tau = 6.0
    dose = 1.0
    cl = 0.25
    v1 = 0.5
    q = 0.5
    v2 = 1.0
    coefficient_a, alpha, coefficient_b, beta = two_compartment_macro_constants(
        dose_mg=dose,
        cl_l_h=cl,
        v1_l=v1,
        q_l_h=q,
        v2_l=v2,
    )
    cmax = coefficient_a / (1.0 - math.exp(-alpha * tau)) + coefficient_b / (1.0 - math.exp(-beta * tau))
    result = simulate_exposure(
        ExposureSimulationRequest.model_validate(
            {
                "simulation_id": "SIM-2C-SS",
                "context_of_use": "Research check of a steady-state two-compartment IV bolus interval.",
                "model_id": "two_compartment_iv_bolus",
                "dose": dose,
                "dose_unit": "mg",
                "clearance": cl,
                "clearance_unit": "L/h",
                "central_volume": v1,
                "central_volume_unit": "L",
                "intercompartmental_clearance": q,
                "intercompartmental_clearance_unit": "L/h",
                "peripheral_volume": v2,
                "peripheral_volume_unit": "L",
                "dosing_interval": tau,
                "dosing_interval_unit": "h",
                "steady_state": True,
                "times_h": [0, tau],
            }
        )
    )
    assert result.cmax == pytest.approx(cmax, rel=1e-12)
    assert result.auc_tau == pytest.approx(dose / cl)
    assert result.auc_inf is None
    assert result.tmax_h == 0.0


def test_regimen_scope_failures() -> None:
    with pytest.raises(ValidationError, match="dosing_interval"):
        _iv_request(n_doses=2)
    with pytest.raises(ValidationError, match="dosing_interval"):
        _iv_request(steady_state=True)
    with pytest.raises(ValidationError, match="n_doses must stay 1"):
        _iv_request(steady_state=True, n_doses=3, dosing_interval=8.0, dosing_interval_unit="h")
    with pytest.raises(AnalysisError, match="shorter than the dosing interval"):
        simulate_exposure(
            ExposureSimulationRequest.model_validate(
                {
                    "simulation_id": "SIM-OVERLAP",
                    "context_of_use": "Research request whose infusion lasts longer than the dosing interval.",
                    "model_id": "one_compartment_iv_infusion",
                    "dose": 1.0,
                    "dose_unit": "mg",
                    "clearance": 0.2,
                    "clearance_unit": "L/h",
                    "volume": 1.0,
                    "volume_unit": "L",
                    "infusion_duration": 2.0,
                    "infusion_duration_unit": "h",
                    "dosing_interval": 1.0,
                    "dosing_interval_unit": "h",
                    "n_doses": 2,
                    "times_h": [0, 1],
                }
            )
        )
    with pytest.raises(AnalysisError, match="one dosing interval"):
        simulate_exposure(
            _iv_request(
                steady_state=True,
                dosing_interval=2.0,
                dosing_interval_unit="h",
                times_h=[0, 3],
                threshold_concentration=None,
                threshold_unit=None,
            )
        )
