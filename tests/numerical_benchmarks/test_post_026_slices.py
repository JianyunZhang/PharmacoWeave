"""Closed forms that follow the weighted two-stage summary."""

import json
import math
from pathlib import Path

import pytest

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pd_models.effect_compartment import simulate_effect_compartment
from pharmacoweave.pk_models.analytical import (
    three_compartment_iv_bolus_concentration,
    three_compartment_macro_constants,
)
from pharmacoweave.pk_models.fit_three_compartment import (
    fit_three_compartment_iv_bolus_profile,
    recover_three_compartment_iv_bolus,
)
from pharmacoweave.qsp.inhibited_turnover import simulate_constant_inhibition_turnover
from pharmacoweave.qsp.turnover import simulate_linear_turnover

_REPO = Path(__file__).resolve().parents[2]
_INHIBITED = json.loads((_REPO / "tests" / "regression" / "inhibited_turnover.json").read_text(encoding="utf-8"))
_EFFECT = json.loads((_REPO / "tests" / "regression" / "effect_compartment.json").read_text(encoding="utf-8"))
_THREE = json.loads((_REPO / "tests" / "regression" / "three_compartment_iv_bolus.json").read_text(encoding="utf-8"))


def test_constant_inhibition_matches_the_hand_formula_and_the_uninhibited_curve() -> None:
    inhibited = simulate_constant_inhibition_turnover(
        kin=_INHIBITED["kin"],
        kout=_INHIBITED["kout"],
        r0=_INHIBITED["r0"],
        inhibition=_INHIBITED["inhibition"],
        times_h=[0.0, math.log(2.0)],
        context_of_use="Research check of constant inhibition.",
    )
    assert inhibited.evidence_grade == "predicted"
    assert inhibited.steady_state == pytest.approx(_INHIBITED["steady_state"], rel=1e-12, abs=1e-12)
    assert inhibited.curve[0].response == pytest.approx(_INHIBITED["response_at_0"], abs=1e-12)
    assert inhibited.curve[1].response == pytest.approx(_INHIBITED["response_at_ln2"], rel=1e-12, abs=1e-12)
    plain = simulate_linear_turnover(
        kin=_INHIBITED["kin"],
        kout=_INHIBITED["kout"],
        r0=_INHIBITED["r0"],
        times_h=[0.0, 1.0, 2.0],
        context_of_use="Research check of constant inhibition.",
    )
    zero = simulate_constant_inhibition_turnover(
        kin=_INHIBITED["kin"],
        kout=_INHIBITED["kout"],
        r0=_INHIBITED["r0"],
        inhibition=0.0,
        times_h=[0.0, 1.0, 2.0],
        context_of_use="Research check of constant inhibition.",
    )
    for left, right in zip(zero.curve, plain.curve, strict=True):
        assert left.response == pytest.approx(right.response, rel=1e-12, abs=1e-12)
    assert zero.steady_state == pytest.approx(plain.steady_state, rel=1e-12, abs=1e-12)
    with pytest.raises(AnalysisError, match="Inhibition"):
        simulate_constant_inhibition_turnover(
            kin=2, kout=1, r0=0, inhibition=1, times_h=[0, 1], context_of_use="research"
        )
    with pytest.raises(AnalysisError, match="does not estimate"):
        simulate_constant_inhibition_turnover(
            kin=2, kout=1, r0=0, inhibition=0.5, times_h=[0, 1], context_of_use="research", estimate_parameters=True
        )


def test_effect_compartment_matches_the_hand_formula() -> None:
    curve = simulate_effect_compartment(
        dose_mg=_EFFECT["dose_mg"],
        volume_l=_EFFECT["volume_l"],
        clearance_l_h=_EFFECT["clearance_l_h"],
        ke0_per_h=_EFFECT["ke0_per_h"],
        times_h=[0.0, math.log(2.0)],
        context_of_use="Research check of an effect compartment.",
    )
    assert curve.evidence_grade == "predicted"
    assert curve.curve[0].effect_site_concentration_mg_l == pytest.approx(_EFFECT["effect_site_at_0"], abs=1e-12)
    assert curve.curve[0].plasma_concentration_mg_l == pytest.approx(1.0, rel=1e-12, abs=1e-12)
    assert curve.curve[1].effect_site_concentration_mg_l == pytest.approx(
        _EFFECT["effect_site_at_ln2_mg_l"], rel=1e-12, abs=1e-12
    )
    equal = simulate_effect_compartment(
        dose_mg=1,
        volume_l=1,
        clearance_l_h=1,
        ke0_per_h=_EFFECT["equal_rate_ke0_per_h"],
        times_h=[_EFFECT["equal_rate_time_h"]],
        context_of_use="Research check of an effect compartment.",
    )
    assert equal.curve[0].effect_site_concentration_mg_l == pytest.approx(math.exp(-1.0), rel=1e-12, abs=1e-12)
    assert equal.constitutes_dosing_recommendation is False
    with pytest.raises(AnalysisError, match="ke0"):
        simulate_effect_compartment(
            dose_mg=1, volume_l=1, clearance_l_h=1, ke0_per_h=0, times_h=[0, 1], context_of_use="research"
        )
    with pytest.raises(AnalysisError, match="does not estimate"):
        simulate_effect_compartment(
            dose_mg=1,
            volume_l=1,
            clearance_l_h=1,
            ke0_per_h=2,
            times_h=[0, 1],
            context_of_use="research",
            estimate_parameters=True,
        )


def test_three_compartment_closed_form_and_noise_free_recovery() -> None:
    dose = _THREE["dose_mg"]
    kwargs = {
        "dose": dose,
        "dose_unit": "mg",
        "clearance": _THREE["clearance_l_h"],
        "clearance_unit": "L/h",
        "central_volume": _THREE["central_volume_l"],
        "central_volume_unit": "L",
        "q2": _THREE["q2_l_h"],
        "q2_unit": "L/h",
        "peripheral_volume_2": _THREE["v2_l"],
        "peripheral_volume_2_unit": "L",
        "q3": _THREE["q3_l_h"],
        "q3_unit": "L/h",
        "peripheral_volume_3": _THREE["v3_l"],
        "peripheral_volume_3_unit": "L",
        "time_unit": "h",
        "concentration_unit": "mg/L",
    }
    assert three_compartment_iv_bolus_concentration(time=0, **kwargs) == pytest.approx(_THREE["c0_mg_l"], rel=1e-12)
    macro = three_compartment_macro_constants(
        dose_mg=dose,
        cl_l_h=_THREE["clearance_l_h"],
        v1_l=_THREE["central_volume_l"],
        q2_l_h=_THREE["q2_l_h"],
        v2_l=_THREE["v2_l"],
        q3_l_h=_THREE["q3_l_h"],
        v3_l=_THREE["v3_l"],
    )
    coefficient_auc = macro[0] / macro[1] + macro[2] / macro[3] + macro[4] / macro[5]
    assert coefficient_auc == pytest.approx(_THREE["auc_inf_mg_h_l"], rel=1e-12, abs=1e-12)
    times = [0, 0.05, 0.1, 0.2, 0.4, 0.8, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 48, 72]
    concentrations = [three_compartment_iv_bolus_concentration(time=time, **kwargs) for time in times]
    recovered = recover_three_compartment_iv_bolus(dose_mg=dose, times_h=times, concentrations_mg_l=concentrations)
    tolerance = _THREE["recovery_relative_error"]
    assert recovered.converged is True
    assert recovered.clearance_l_h == pytest.approx(_THREE["clearance_l_h"], rel=tolerance)
    assert recovered.central_volume_l == pytest.approx(_THREE["central_volume_l"], rel=tolerance)
    assert recovered.q2_l_h == pytest.approx(_THREE["q2_l_h"], rel=tolerance)
    assert recovered.v2_l == pytest.approx(_THREE["v2_l"], rel=tolerance)
    assert recovered.q3_l_h == pytest.approx(_THREE["q3_l_h"], rel=tolerance)
    assert recovered.v3_l == pytest.approx(_THREE["v3_l"], rel=tolerance)
    failed = recover_three_compartment_iv_bolus(dose_mg=dose, times_h=[0, 1, 2], concentrations_mg_l=[1, 0.5, 0.2])
    assert failed.converged is False
    assert failed.clearance_l_h is None


def test_three_compartment_profile_fit_reports_intervals_and_fails_closed() -> None:
    times = [0, 0.05, 0.1, 0.2, 0.4, 0.8, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 48, 72]
    points = tuple(
        QuantifiedPoint(
            observation_id=f"OBS-{index:02d}",
            time_h=time,
            concentration_mg_l=three_compartment_iv_bolus_concentration(
                dose=_THREE["dose_mg"],
                dose_unit="mg",
                clearance=_THREE["clearance_l_h"],
                clearance_unit="L/h",
                central_volume=_THREE["central_volume_l"],
                central_volume_unit="L",
                q2=_THREE["q2_l_h"],
                q2_unit="L/h",
                peripheral_volume_2=_THREE["v2_l"],
                peripheral_volume_2_unit="L",
                q3=_THREE["q3_l_h"],
                q3_unit="L/h",
                peripheral_volume_3=_THREE["v3_l"],
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
        dose_mg=_THREE["dose_mg"],
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    fitted = fit_three_compartment_iv_bolus_profile(profile, weight="uniform")
    assert fitted.diagnostic.status is DiagnosticStatus.PASS
    by_name = {estimate.parameter: estimate for estimate in fitted.estimates}
    tolerance = _THREE["recovery_relative_error"]
    assert by_name["CL"].value == pytest.approx(_THREE["clearance_l_h"], rel=tolerance)
    assert by_name["V1"].value == pytest.approx(_THREE["central_volume_l"], rel=tolerance)
    assert by_name["Q3"].value == pytest.approx(_THREE["q3_l_h"], rel=tolerance)
    assert by_name["CL"].standard_error is not None
    assert by_name["CL"].ci_lower is not None and by_name["CL"].ci_upper is not None
    short = ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="DRUG-01",
        dose_event_id="DOSE-01",
        route=Route.IV_BOLUS,
        dose_mg=_THREE["dose_mg"],
        points=points[:3],
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    failed = fit_three_compartment_iv_bolus_profile(short, weight="uniform")
    assert failed.diagnostic.status is DiagnosticStatus.FAIL
    assert failed.estimates == ()
