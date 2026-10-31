"""Numerical benchmarks for 1-compartment extravascular, 2-compartment IV bolus, and molar PK."""

import pytest
from examples.iv_bolus_case import build_study as build_iv_one_comp_study
from examples.oral_absorption_case import build_noisy_oral_profile, build_oral_study
from examples.two_compartment_iv_case import (
    DOSE_MG,
    TRUE_CL_L_H,
    TRUE_Q_L_H,
    TRUE_V1_L,
    TRUE_V2_L,
    build_noisy_two_compartment_profile,
    build_two_compartment_study,
)

from pharmacoweave.core.entities import ObservationStatus, Route
from pharmacoweave.core.modeling import DiagnosticStatus, ParameterSemantic
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.errors import AnalysisError, ValidationFailed
from pharmacoweave.pk_models.analytical import one_compartment_oral_concentration
from pharmacoweave.pk_models.fit_extravascular import (
    fit_one_compartment_extravascular,
    fit_one_compartment_extravascular_profile,
)
from pharmacoweave.pk_models.fit_two_compartment import (
    fit_two_compartment_iv_bolus,
    fit_two_compartment_iv_bolus_profile,
)


def test_quiet_oral_fit_recovers_ka_cl_over_f_and_v_over_f() -> None:
    fits = fit_one_compartment_extravascular(build_oral_study())
    assert len(fits) == 2
    by_subject = {item.subject_id: item for item in fits}
    expected = {
        "SUBJ-01": {"ka": 1.50, "CL/F": 0.0625, "V/F": 0.25},
        "SUBJ-02": {"ka": 1.50, "CL/F": 0.0750, "V/F": 0.30},
    }
    for subject_id, result in by_subject.items():
        assert result.fit.converged is True
        assert result.fit.diagnostic_status is DiagnosticStatus.PASS
        assert result.fit.limitations
        estimates = {est.parameter: est for est in result.estimates}
        assert set(estimates) == {"ka", "CL/F", "V/F"}
        assert estimates["ka"].semantic is ParameterSemantic.KA
        assert estimates["CL/F"].semantic is ParameterSemantic.CL_OVER_F
        assert estimates["V/F"].semantic is ParameterSemantic.V_OVER_F
        for name, target in expected[subject_id].items():
            est = estimates[name]
            assert est.value == pytest.approx(target, rel=1e-5)
            assert est.standard_error is not None and est.standard_error == pytest.approx(0.0, abs=1e-4)
            assert est.ci_lower is not None and est.ci_upper is not None
            assert 0 < est.ci_lower <= est.value <= est.ci_upper


def test_noisy_oral_fit_confidence_intervals_cover_true_parameters() -> None:
    profile = build_noisy_oral_profile()
    result = fit_one_compartment_extravascular_profile(profile, weight="1/y^2")
    assert result.fit.converged is True
    assert result.fit.diagnostic_status is DiagnosticStatus.WARNING
    estimates = {est.parameter: est for est in result.estimates}
    targets = {"ka": 1.50, "CL/F": 0.25, "V/F": 1.0}
    for name, target in targets.items():
        est = estimates[name]
        assert est.value == pytest.approx(target, rel=0.12)
        assert est.standard_error is not None and est.standard_error > 0
        assert est.ci_lower is not None and est.ci_upper is not None
        assert est.ci_lower <= target <= est.ci_upper


def test_oral_fit_fails_without_absorption_phase_points() -> None:
    declining_only = ConcentrationProfile(
        subject_id="SUBJ-NO-ABS",
        analyte_id="PW-0001-PARENT",
        dose_event_id="DOSE-1",
        route=Route.ORAL,
        dose_mg=2.0,
        points=(
            QuantifiedPoint("OBS-1", 1.5, 1.40),
            QuantifiedPoint("OBS-2", 2.0, 1.33),
            QuantifiedPoint("OBS-3", 4.0, 0.88),
            QuantifiedPoint("OBS-4", 8.0, 0.32),
            QuantifiedPoint("OBS-5", 12.0, 0.12),
        ),
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    result = fit_one_compartment_extravascular_profile(declining_only)
    assert result.fit.diagnostic_status is DiagnosticStatus.FAIL
    assert result.fit.converged is False
    assert result.estimates == ()
    assert result.predicted_mg_l == ()


def test_oral_fit_fails_when_ka_and_ke_confluence_is_non_identifiable() -> None:
    times_h = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
    points = tuple(
        QuantifiedPoint(
            observation_id=f"OBS-CONF-{idx}",
            time_h=time_h,
            concentration_mg_l=one_compartment_oral_concentration(
                dose=2.0,
                dose_unit="mg",
                ka=0.40,
                ka_unit="1/h",
                apparent_volume=1.0,
                volume_unit="L",
                apparent_clearance=0.40,
                clearance_unit="L/h",
                time=time_h,
                time_unit="h",
                concentration_unit="mg/L",
            ),
        )
        for idx, time_h in enumerate(times_h, start=1)
    )
    confluent = ConcentrationProfile(
        subject_id="SUBJ-CONFLUENT",
        analyte_id="PW-0001-PARENT",
        dose_event_id="DOSE-CONF",
        route=Route.ORAL,
        dose_mg=2.0,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
    )
    result = fit_one_compartment_extravascular_profile(confluent)
    assert result.fit.diagnostic_status is DiagnosticStatus.FAIL
    assert result.fit.converged is False
    assert result.estimates == ()
    assert result.predicted_mg_l == ()


def test_quiet_two_compartment_fit_recovers_cl_v1_q_and_v2() -> None:
    fits = fit_two_compartment_iv_bolus(build_two_compartment_study())
    assert len(fits) == 1
    result = fits[0]
    assert result.fit.converged is True
    assert result.fit.diagnostic_status is DiagnosticStatus.PASS
    assert result.fit.limitations
    estimates = {est.parameter: est for est in result.estimates}
    assert set(estimates) == {"CL", "V1", "Q", "V2"}
    assert estimates["CL"].semantic is ParameterSemantic.CL
    assert estimates["V1"].semantic is ParameterSemantic.V1
    assert estimates["Q"].semantic is ParameterSemantic.Q
    assert estimates["V2"].semantic is ParameterSemantic.V2
    targets = {"CL": TRUE_CL_L_H, "V1": TRUE_V1_L, "Q": TRUE_Q_L_H, "V2": TRUE_V2_L}
    for name, target in targets.items():
        est = estimates[name]
        assert est.value == pytest.approx(target, rel=1e-5)
        assert est.standard_error is not None and est.standard_error == pytest.approx(0.0, abs=1e-4)
        assert est.ci_lower is not None and est.ci_upper is not None
        assert 0 < est.ci_lower <= est.value <= est.ci_upper


def test_noisy_two_compartment_fit_confidence_intervals_cover_true_parameters() -> None:
    profile = build_noisy_two_compartment_profile()
    result = fit_two_compartment_iv_bolus_profile(profile, weight="1/y^2")
    assert result.fit.converged is True
    assert result.fit.diagnostic_status is DiagnosticStatus.WARNING
    estimates = {est.parameter: est for est in result.estimates}
    targets = {"CL": TRUE_CL_L_H, "V1": TRUE_V1_L, "Q": TRUE_Q_L_H, "V2": TRUE_V2_L}
    for name, target in targets.items():
        est = estimates[name]
        assert est.value == pytest.approx(target, rel=0.15)
        assert est.standard_error is not None and est.standard_error > 0
        assert est.ci_lower is not None and est.ci_upper is not None
        assert est.ci_lower <= target <= est.ci_upper
    assert DOSE_MG == 10.0


def test_two_compartment_fit_rejects_monoexponential_curve() -> None:
    fits = fit_two_compartment_iv_bolus(build_iv_one_comp_study())
    assert len(fits) == 2
    for result in fits:
        assert result.fit.diagnostic_status is DiagnosticStatus.FAIL
        assert result.fit.converged is False
        assert result.estimates == ()
        assert result.predicted_mg_l == ()


def test_molar_study_converts_via_explicit_molecular_weight_and_fails_without_it() -> None:
    mass_study = build_oral_study()
    mw = 250.0  # g/mol: 1 mg = 4 umol; 1 ug/mL = 4 umol/L = 4000 nmol/L
    molar_doses = [
        dose.model_copy(
            update={
                "dose_value": dose.dose_value * 1000.0 / mw,
                "dose_unit": "umol",
                "body_weight_normalized_value": None,
                "body_weight_normalized_unit": None,
            }
        )
        for dose in mass_study.dose_events
    ]
    molar_assays = [
        assay.model_copy(
            update={
                "lloq": assay.lloq * 1000.0 / mw,
                "uloq": assay.uloq * 1000.0 / mw,
                "limit_unit": "umol/L",
            }
        )
        for assay in mass_study.assays
    ]
    molar_observations = []
    for obs in mass_study.observations:
        if obs.status is ObservationStatus.QUANTIFIED and obs.numeric_value is not None:
            nmol_l = obs.numeric_value * 1_000_000.0 / mw
            text = f"{nmol_l:.6f}"
            molar_observations.append(
                obs.model_copy(
                    update={
                        "unit": "nmol/L",
                        "numeric_value": nmol_l,
                        "raw_value": text,
                        "reported_value": text,
                    }
                )
            )
        else:
            molar_observations.append(obs.model_copy(update={"unit": "nmol/L"}))

    molar_study = mass_study.model_copy(
        update={
            "dose_events": molar_doses,
            "assays": molar_assays,
            "observations": molar_observations,
        }
    )
    mass_nca = {row.subject_id: row for row in mass_study.run_nca("configs/nca_single_dose_oral.yaml")}
    molar_nca = {row.subject_id: row for row in molar_study.run_nca("configs/nca_single_dose_oral.yaml")}
    for subject_id, expected in mass_nca.items():
        actual = molar_nca[subject_id]
        assert actual.auc_last == pytest.approx(expected.auc_last, rel=1e-10)
        assert actual.auc_inf == pytest.approx(expected.auc_inf, rel=1e-10)
        assert actual.clearance == pytest.approx(expected.clearance, rel=1e-10)
        assert actual.volume == pytest.approx(expected.volume, rel=1e-10)

    molar_fits = {item.subject_id: item for item in fit_one_compartment_extravascular(molar_study)}
    assert molar_fits["SUBJ-01"].fit.diagnostic_status is DiagnosticStatus.PASS
    estimates = {est.parameter: est.value for est in molar_fits["SUBJ-01"].estimates}
    assert estimates["ka"] == pytest.approx(1.50, rel=1e-4)
    assert estimates["CL/F"] == pytest.approx(0.0625, rel=1e-4)
    assert estimates["V/F"] == pytest.approx(0.25, rel=1e-4)

    missing_mw_compounds = [
        comp.model_copy(update={"molecular_weight": None, "molecular_weight_unit": None})
        for comp in molar_study.compounds
    ]
    missing_mw_study = molar_study.model_copy(update={"compounds": missing_mw_compounds})
    with pytest.raises((AnalysisError, ValidationFailed)):
        fit_one_compartment_extravascular(missing_mw_study)
