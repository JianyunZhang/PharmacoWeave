"""Closed-form checks for exploratory single-dose exposure simulation."""

import math

import pytest
from pydantic import ValidationError

from pharmacoweave.errors import ContractViolation, PlannedCapabilityError
from pharmacoweave.pk_models.analytical import bateman_mg_l, iv_infusion_concentration
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
