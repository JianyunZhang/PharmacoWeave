"""Weighted nonlinear least squares for the one-compartment IV bolus."""

import math

import pytest
from examples.iv_bolus_noisy import (
    EXAMPLE_DIR,
    RECOVERY_RELATIVE_TOLERANCE,
    TRUE_CLEARANCE_L_H,
    TRUE_VOLUME_L,
    build_study,
)

from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.fit_iv import fit_iv_bolus
from pharmacoweave.pk_models.wnls import fit_iv_bolus_profile_wnls, fit_iv_bolus_wnls


def _profile(
    times: list[float],
    concentrations: list[float] | None = None,
    *,
    dose_mg: float = 0.25,
) -> ConcentrationProfile:
    if concentrations is None:
        concentrations = [math.exp(-0.2 * time_h) for time_h in times]
    points = tuple(
        QuantifiedPoint(observation_id=f"OBS-{index}", time_h=time_h, concentration_mg_l=concentration)
        for index, (time_h, concentration) in enumerate(zip(times, concentrations, strict=True))
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


def _values(fit) -> dict[str, float]:
    return {estimate.parameter: estimate.value for estimate in fit.estimates}


@pytest.mark.parametrize("weight", ["1/y", "1/y^2"])
def test_noise_free_curve_recovers_clearance_and_volume(weight: str) -> None:
    fit = fit_iv_bolus_profile_wnls(_profile([0.5, 1, 2, 4, 8, 12]), weight=weight)
    assert fit.fit.converged is True
    assert fit.fit.diagnostic_status is DiagnosticStatus.PASS
    values = _values(fit)
    assert values["CL"] == pytest.approx(0.05, rel=1e-8)
    assert values["V"] == pytest.approx(0.25, rel=1e-8)
    assert fit.fit.limitations


@pytest.mark.parametrize("weight", ["1/y", "1/y^2"])
def test_noisy_example_recovers_clearance_and_volume(weight: str) -> None:
    fits = fit_iv_bolus_wnls(build_study(), weight=weight)
    assert len(fits) == 1
    fit = fits[0]
    assert fit.fit.converged is True
    assert fit.diagnostic.status is DiagnosticStatus.WARNING
    assert fit.estimates
    values = _values(fit)
    assert values["CL"] == pytest.approx(TRUE_CLEARANCE_L_H, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert values["V"] == pytest.approx(TRUE_VOLUME_L, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert all("T1440" not in observation_id and "PRE" not in observation_id for observation_id in fit.observation_ids)
    assert min(fit.observed_mg_l) > 0


def test_checked_in_noisy_study_validates_and_recovers() -> None:
    study = PharmacologyStudy.from_files(EXAMPLE_DIR, mode="strict")
    fit = fit_iv_bolus_wnls(study, weight="1/y^2")[0]
    values = _values(fit)
    assert values["CL"] == pytest.approx(TRUE_CLEARANCE_L_H, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert values["V"] == pytest.approx(TRUE_VOLUME_L, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert study.observations[-1].status.value == "BLQ"
    assert study.observations[-1].numeric_value is None


def test_log_linear_path_still_fits_the_noisy_example() -> None:
    fit = fit_iv_bolus(build_study())[0]
    assert fit.fit.estimation_method == "log_linear_ols"
    assert fit.fit.converged is True
    values = _values(fit)
    assert values["CL"] == pytest.approx(TRUE_CLEARANCE_L_H, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert values["V"] == pytest.approx(TRUE_VOLUME_L, rel=RECOVERY_RELATIVE_TOLERANCE)


def test_increasing_concentrations_do_not_return_clearance() -> None:
    profile = _profile([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], dose_mg=1.0)
    fit = fit_iv_bolus_profile_wnls(profile, weight="1/y^2")
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL
    assert fit.fit.converged is False
    assert fit.predicted_mg_l == ()
    assert not any(estimate.parameter == "CL" for estimate in fit.estimates)


def test_repeated_times_are_not_identifiable() -> None:
    profile = _profile([1.0, 1.0, 1.0], [1.0, 0.8, 0.5], dose_mg=1.0)
    fit = fit_iv_bolus_profile_wnls(profile, weight="1/y")
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL
    assert any("not identifiable" in item for item in fit.fit.limitations)


def test_log_linear_repeated_times_do_not_raise() -> None:
    from pharmacoweave.pk_models.fit_iv import fit_iv_bolus_profile

    fit = fit_iv_bolus_profile(_profile([2.0, 2.0, 2.0], [1.0, 0.5, 0.25], dose_mg=1.0))
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL


def test_two_points_do_not_return_clearance() -> None:
    fit = fit_iv_bolus_profile_wnls(_profile([1.0, 2.0]), weight="1/y^2")
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL


def test_unknown_weight_is_rejected() -> None:
    with pytest.raises(AnalysisError, match="weight"):
        fit_iv_bolus_profile_wnls(_profile([0.5, 1, 2, 4]), weight="1/sqrt(y)")
