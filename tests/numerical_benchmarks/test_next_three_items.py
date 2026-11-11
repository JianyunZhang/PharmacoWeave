"""Concentration-driven turnover, single-dose three-compartment exposure, and audit fit choice."""

import json
import math
from pathlib import Path

import pytest
from pydantic import ValidationError

from pharmacoweave.cli import main
from pharmacoweave.core.entities import Route
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.errors import AnalysisError
from pharmacoweave.integrations.pharmacology_evidence import evidence_from_pk_fit
from pharmacoweave.pk_models.analytical import three_compartment_iv_bolus_concentration
from pharmacoweave.pk_models.fit_three_compartment import fit_three_compartment_iv_bolus_profile
from pharmacoweave.qsp.concentration_turnover import simulate_concentration_driven_turnover
from pharmacoweave.simulation.exposure import ExposureSimulationRequest, simulate_exposure

_REPO = Path(__file__).resolve().parents[2]
_TURN = json.loads((_REPO / "tests" / "regression" / "concentration_turnover.json").read_text(encoding="utf-8"))
_EXPOSURE = json.loads((_REPO / "tests" / "regression" / "three_compartment_exposure.json").read_text(encoding="utf-8"))


def test_concentration_driven_turnover_matches_the_hand_formula() -> None:
    equal = simulate_concentration_driven_turnover(
        kin=_TURN["kin"],
        kout=_TURN["equal_kout"],
        r0=0,
        dose_mg=_TURN["dose_mg"],
        volume_l=_TURN["volume_l"],
        clearance_l_h=_TURN["clearance_l_h"],
        times_h=[0, _TURN["equal_time_h"]],
        context_of_use="Research check of concentration-driven production.",
    )
    assert equal.evidence_grade == "predicted"
    assert equal.steady_state == 0
    assert equal.curve[0].response == pytest.approx(0, abs=1e-12)
    assert equal.curve[0].plasma_concentration_mg_l == pytest.approx(1, rel=1e-12, abs=1e-12)
    assert equal.curve[1].response == pytest.approx(math.exp(-1.0), rel=1e-12, abs=1e-12)
    distinct = simulate_concentration_driven_turnover(
        kin=_TURN["distinct_kin"],
        kout=_TURN["distinct_kout"],
        r0=1,
        dose_mg=_TURN["dose_mg"],
        volume_l=_TURN["volume_l"],
        clearance_l_h=_TURN["clearance_l_h"],
        times_h=[0, math.log(2.0)],
        context_of_use="Research check of concentration-driven production.",
    )
    assert distinct.curve[0].response == pytest.approx(1, rel=1e-12, abs=1e-12)
    assert distinct.curve[1].response == pytest.approx(
        math.exp(-2.0 * math.log(2.0)) + _TURN["distinct_response_at_ln2"], rel=1e-12, abs=1e-12
    )
    with pytest.raises(AnalysisError, match="does not estimate"):
        simulate_concentration_driven_turnover(
            kin=1,
            kout=1,
            r0=0,
            dose_mg=1,
            volume_l=1,
            clearance_l_h=1,
            times_h=[0, 1],
            context_of_use="research",
            estimate_parameters=True,
        )


def test_three_compartment_single_dose_exposure_and_repeated_rejection() -> None:
    result = simulate_exposure(
        ExposureSimulationRequest(
            simulation_id="SIM-3C",
            context_of_use="Research single-dose three-compartment exposure check.",
            model_id="three_compartment_iv_bolus",
            dose=_EXPOSURE["dose_mg"],
            dose_unit="mg",
            clearance=_EXPOSURE["clearance_l_h"],
            clearance_unit="L/h",
            central_volume=_EXPOSURE["central_volume_l"],
            central_volume_unit="L",
            intercompartmental_clearance=_EXPOSURE["q2_l_h"],
            intercompartmental_clearance_unit="L/h",
            peripheral_volume=_EXPOSURE["v2_l"],
            peripheral_volume_unit="L",
            q3=_EXPOSURE["q3_l_h"],
            q3_unit="L/h",
            peripheral_volume_3=_EXPOSURE["v3_l"],
            peripheral_volume_3_unit="L",
            times_h=[0, 1, 4],
        )
    )
    assert result.evidence_type == "predicted"
    assert result.regimen == "single_dose"
    assert result.terminal_rate_name == "gamma"
    assert result.curve[0].concentration == pytest.approx(_EXPOSURE["c0_mg_l"], rel=1e-12, abs=1e-12)
    assert result.auc_inf == pytest.approx(_EXPOSURE["auc_inf_mg_h_l"], rel=1e-12, abs=1e-12)
    direct = three_compartment_iv_bolus_concentration(
        dose=_EXPOSURE["dose_mg"],
        dose_unit="mg",
        clearance=_EXPOSURE["clearance_l_h"],
        clearance_unit="L/h",
        central_volume=_EXPOSURE["central_volume_l"],
        central_volume_unit="L",
        q2=_EXPOSURE["q2_l_h"],
        q2_unit="L/h",
        peripheral_volume_2=_EXPOSURE["v2_l"],
        peripheral_volume_2_unit="L",
        q3=_EXPOSURE["q3_l_h"],
        q3_unit="L/h",
        peripheral_volume_3=_EXPOSURE["v3_l"],
        peripheral_volume_3_unit="L",
        time=1,
        time_unit="h",
        concentration_unit="mg/L",
    )
    assert result.curve[1].concentration == pytest.approx(direct, rel=1e-12, abs=1e-12)
    with pytest.raises(ValidationError, match="single-dose only"):
        ExposureSimulationRequest(
            simulation_id="SIM-3C-REPEAT",
            context_of_use="Research single-dose three-compartment exposure check.",
            model_id="three_compartment_iv_bolus",
            dose=1,
            dose_unit="mg",
            clearance=0.5,
            clearance_unit="L/h",
            central_volume=1,
            central_volume_unit="L",
            intercompartmental_clearance=1.5,
            intercompartmental_clearance_unit="L/h",
            peripheral_volume=4,
            peripheral_volume_unit="L",
            q3=0.25,
            q3_unit="L/h",
            peripheral_volume_3=10,
            peripheral_volume_3_unit="L",
            times_h=[0, 1],
            n_doses=2,
            dosing_interval=12,
            dosing_interval_unit="h",
        )


def test_audit_fit_model_accepts_three_compartment_and_exports_clearance(capsys: pytest.CaptureFixture[str]) -> None:
    times = [0, 0.05, 0.1, 0.2, 0.4, 0.8, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 48, 72]
    points = tuple(
        QuantifiedPoint(
            observation_id=f"OBS-{index:02d}",
            time_h=time,
            concentration_mg_l=three_compartment_iv_bolus_concentration(
                dose=1,
                dose_unit="mg",
                clearance=0.5,
                clearance_unit="L/h",
                central_volume=1,
                central_volume_unit="L",
                q2=1.5,
                q2_unit="L/h",
                peripheral_volume_2=4,
                peripheral_volume_2_unit="L",
                q3=0.25,
                q3_unit="L/h",
                peripheral_volume_3=10,
                peripheral_volume_3_unit="L",
                time=time,
                time_unit="h",
                concentration_unit="mg/L",
            ),
        )
        for index, time in enumerate(times, start=1)
    )
    profile = ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="DRUG-01",
        dose_event_id="DOSE-01",
        route=Route.IV_BOLUS,
        dose_mg=1,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    fitted = fit_three_compartment_iv_bolus_profile(profile, weight="uniform")
    evidence = evidence_from_pk_fit(fitted, run_id="audit-3c")
    clearance = next(item for item in evidence if item.claim_type == "CL")
    assert clearance.numeric_value == pytest.approx(0.5, rel=1e-4)
    with pytest.raises(SystemExit) as caught:
        main(["audit-predictions", "--help"])
    assert caught.value.code in {0, None}
    assert "three_compartment_iv_bolus" in capsys.readouterr().out
