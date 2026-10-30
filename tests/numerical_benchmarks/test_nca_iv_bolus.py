"""Numerical checks for single-dose IV bolus NCA and the log-linear fit."""

import math

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import auc_inf_is_reportable
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.fit_iv import fit_iv_bolus, fit_iv_bolus_profile
from pharmacoweave.pk_nca.calculate import nca_iv_bolus
from pharmacoweave.pk_nca.config import NcaConfig


def _config() -> NcaConfig:
    return NcaConfig.model_validate(
        {
            "analysis": "nca",
            "design_route": "iv_bolus",
            "dosing": "single_dose",
            "sampling": "rich_serial",
            "auc_method": "linear_up_log_down",
            "lambda_z_min_points": 3,
            "lambda_z_min_adjusted_r2": 0.90,
            "auc_inf_max_extrapolated_fraction": 0.20,
            "blq_rule": "keep_as_censored",
        }
    )


def _profile(times: list[float], *, kel: float = 0.2, dose_mg: float = 0.25) -> ConcentrationProfile:
    points = tuple(
        QuantifiedPoint(observation_id=f"OBS-{index}", time_h=time_h, concentration_mg_l=math.exp(-kel * time_h))
        for index, time_h in enumerate(times)
    )
    return ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="AN-1",
        dose_event_id="DOSE-1",
        route=Route.IV_BOLUS,
        dose_mg=dose_mg,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )


def test_exact_monoexponential_recovers_dose_over_clearance() -> None:
    result = nca_iv_bolus(_profile([0.5, 1, 2, 4, 8, 12]), _config())
    assert result.terminal_phase_qualified
    assert result.lambda_z == pytest.approx(0.2)
    assert result.half_life == pytest.approx(math.log(2.0) / 0.2)
    assert result.auc_inf == pytest.approx(5.0)
    assert result.auc_last == pytest.approx(5.0 * (1.0 - math.exp(-0.2 * 12)))
    assert result.clearance == pytest.approx(0.05)
    assert result.volume == pytest.approx(0.25)
    assert result.clearance_semantic is not None
    assert result.clearance_semantic.value == "CL"
    assert result.volume_semantic is not None
    assert result.volume_semantic.value == "Vz"
    assert auc_inf_is_reportable(result)


def test_two_points_do_not_report_half_life() -> None:
    result = nca_iv_bolus(_profile([1.0, 2.0]), _config())
    assert result.terminal_phase_qualified is False
    assert result.half_life is None
    assert result.auc_inf is None
    assert result.auc_last is not None
    assert result.cmax == pytest.approx(math.exp(-0.2))


def test_early_curve_keeps_auc_inf_and_flags_the_extrapolation() -> None:
    result = nca_iv_bolus(_profile([0.25, 0.5, 0.75]), _config())
    assert result.auc_inf == pytest.approx(5.0)
    assert result.auc_inf_extrap_fraction is not None
    assert result.auc_inf_extrap_fraction > 0.20
    assert auc_inf_is_reportable(result) is False
    assert any("auc_inf_is_reportable is false" in warning for warning in result.warnings)


def test_fixture_matches_the_analytical_targets() -> None:
    results = build_study().run_nca("configs/nca_single_dose_iv.yaml")
    by_subject = {row.subject_id: row for row in results}
    assert by_subject["SUBJ-01"].clearance == pytest.approx(0.05, rel=1e-4)
    assert by_subject["SUBJ-02"].clearance == pytest.approx(0.06, rel=1e-4)
    assert by_subject["SUBJ-01"].volume == pytest.approx(0.25, rel=1e-4)
    assert by_subject["SUBJ-02"].volume == pytest.approx(0.30, rel=1e-4)
    for row in results:
        assert row.auc_inf == pytest.approx(5.0, rel=1e-4)
        assert row.auc_last == pytest.approx(5.0 * (1.0 - math.exp(-2.4)), rel=1e-4)
        assert row.half_life == pytest.approx(math.log(2.0) / 0.2, rel=1e-4)
        assert any("not replaced with a number" in warning for warning in row.warnings)


def test_log_linear_fit_recovers_clearance_and_volume() -> None:
    fit = fit_iv_bolus_profile(_profile([0.5, 1, 2, 4, 8]))
    assert fit.fit.converged is True
    assert fit.fit.diagnostic_status.value == "pass"
    values = {estimate.parameter: estimate.value for estimate in fit.estimates}
    assert values["CL"] == pytest.approx(0.05)
    assert values["V"] == pytest.approx(0.25)
    assert fit.fit.limitations


def test_increasing_concentrations_do_not_become_a_clearance() -> None:
    points = tuple(
        QuantifiedPoint(observation_id=f"OBS-{index}", time_h=time_h, concentration_mg_l=concentration)
        for index, (time_h, concentration) in enumerate([(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)])
    )
    profile = ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="AN-1",
        dose_event_id="DOSE-1",
        route=Route.IV_BOLUS,
        dose_mg=1.0,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    fit = fit_iv_bolus_profile(profile)
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status.value == "fail"
    assert fit.fit.converged is False


def test_fixture_fit_recovers_both_subjects() -> None:
    fits = fit_iv_bolus(build_study())
    values = {fit.fit.fit_id: {item.parameter: item.value for item in fit.estimates} for fit in fits}
    assert values["fit-SUBJ-01-PW-0001-PARENT"]["CL"] == pytest.approx(0.05, rel=1e-4)
    assert values["fit-SUBJ-01-PW-0001-PARENT"]["V"] == pytest.approx(0.25, rel=1e-4)
    assert values["fit-SUBJ-02-PW-0001-PARENT"]["CL"] == pytest.approx(0.06, rel=1e-4)
    assert values["fit-SUBJ-02-PW-0001-PARENT"]["V"] == pytest.approx(0.30, rel=1e-4)


def test_oral_doses_are_rejected() -> None:
    study = build_study()
    doses = [dose.model_copy(update={"route": Route.ORAL}) for dose in study.dose_events]
    edited = study.model_copy(update={"dose_events": doses})
    with pytest.raises(AnalysisError, match="iv_bolus"):
        edited.run_nca("configs/nca_single_dose_iv.yaml")


def test_multiple_doses_are_rejected() -> None:
    study = build_study()
    extra = study.dose_events[0].model_copy(update={"dose_event_id": "DOSE-EXTRA"})
    edited = study.model_copy(update={"dose_events": [*study.dose_events, extra]})
    with pytest.raises(AnalysisError, match="single dose"):
        edited.run_nca("configs/nca_single_dose_iv.yaml")
