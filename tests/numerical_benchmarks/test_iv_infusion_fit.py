"""One-compartment constant-rate IV infusion weighted nonlinear least squares.

Noise-free concentrations come from ``iv_infusion_concentration``. The fit must
recover CL and V, and it must withhold clearance when the infusion duration or
the curvature does not identify them.
"""

import math

import pytest
from examples.iv_infusion_case import (
    CLEARANCE_L_PER_H_PER_KG,
    INFUSION_EXAMPLE_DIR,
    SUBJECTS,
    VOLUME_L_PER_KG,
)

from pharmacoweave.cli import main
from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import DiagnosticStatus, ParameterSemantic
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import AnalysisError
from pharmacoweave.integrations.pharmacology_evidence import evidence_from_pk_fit
from pharmacoweave.pk_models.analytical import iv_infusion_concentration
from pharmacoweave.pk_models.fit_infusion import (
    _concentration_and_log_columns,
    fit_one_compartment_iv_infusion,
    fit_one_compartment_iv_infusion_profile,
)

RECOVERY_RELATIVE_TOLERANCE = 1e-8
DOSE_MG = 1.0
DURATION_H = 0.5
TRUE_CLEARANCE_L_H = 0.2
TRUE_VOLUME_L = 1.0


def _concentration(time_h: float, *, clearance: float = TRUE_CLEARANCE_L_H, volume: float = TRUE_VOLUME_L) -> float:
    return iv_infusion_concentration(
        dose=DOSE_MG,
        dose_unit="mg",
        infusion_duration=DURATION_H,
        infusion_duration_unit="h",
        volume=volume,
        volume_unit="L",
        clearance=clearance,
        clearance_unit="L/h",
        time=time_h,
        time_unit="h",
        concentration_unit="mg/L",
    )


def _profile(
    times: list[float],
    concentrations: list[float] | None = None,
    *,
    dose_mg: float = DOSE_MG,
    duration_h: float | None = DURATION_H,
    route: Route = Route.IV_INFUSION,
    regimen: str = "single_dose",
) -> ConcentrationProfile:
    if concentrations is None:
        concentrations = [_concentration(time_h) for time_h in times]
    points = tuple(
        QuantifiedPoint(observation_id=f"OBS-{index}", time_h=time_h, concentration_mg_l=concentration)
        for index, (time_h, concentration) in enumerate(zip(times, concentrations, strict=True))
    )
    return ConcentrationProfile(
        subject_id="SUBJ-01",
        analyte_id="AN-1",
        dose_event_id="DOSE-1",
        route=route,
        dose_mg=dose_mg,
        points=points,
        n_blq_omitted=0,
        n_above_uloq_omitted=0,
        n_predose_omitted=0,
        infusion_duration_h=duration_h,
        interval_regimen=regimen,  # type: ignore[arg-type]
    )


def _values(fit) -> dict[str, float]:
    return {estimate.parameter: estimate.value for estimate in fit.estimates}


@pytest.mark.parametrize("weight", ["uniform", "1/y", "1/y^2"])
@pytest.mark.parametrize(
    "times",
    [
        [0.25, 0.5, 1.0, 2.0, 4.0, 8.0],
        [1.0, 2.0, 4.0, 8.0],
        [0.1, 0.25, 0.4, 1.0],
    ],
)
def test_noise_free_curve_recovers_clearance_and_volume(weight: str, times: list[float]) -> None:
    fit = fit_one_compartment_iv_infusion_profile(_profile(times), weight=weight)
    assert fit.fit.converged is True
    assert fit.fit.diagnostic_status is DiagnosticStatus.PASS
    assert fit.fit.model_id == "one_compartment_iv_infusion"
    values = _values(fit)
    assert values["CL"] == pytest.approx(TRUE_CLEARANCE_L_H, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert values["V"] == pytest.approx(TRUE_VOLUME_L, rel=RECOVERY_RELATIVE_TOLERANCE)
    assert set(values) == {"CL", "V"}
    for estimate in fit.estimates:
        assert estimate.semantic in {ParameterSemantic.CL, ParameterSemantic.V}
        assert estimate.route is Route.IV_INFUSION
        assert estimate.standard_error == pytest.approx(0.0, abs=1e-8)
        assert estimate.ci_level == 0.95
        assert estimate.ci_lower == pytest.approx(estimate.value, rel=5e-8)
        assert estimate.ci_upper == pytest.approx(estimate.value, rel=5e-8)
    for time_h, predicted in zip(fit.times_h, fit.predicted_mg_l, strict=True):
        reference = iv_infusion_concentration(
            dose=DOSE_MG,
            dose_unit="mg",
            infusion_duration=DURATION_H,
            infusion_duration_unit="h",
            volume=values["V"],
            volume_unit="L",
            clearance=values["CL"],
            clearance_unit="L/h",
            time=time_h,
            time_unit="h",
            concentration_unit="mg/L",
        )
        assert predicted == pytest.approx(reference, rel=1e-12, abs=1e-12)
        assert predicted == pytest.approx(_concentration(time_h), rel=RECOVERY_RELATIVE_TOLERANCE)
    assert any("not a dosing recommendation" in item for item in fit.fit.limitations)
    assert all(estimate.parameter in {"CL", "V"} for estimate in fit.estimates)


def test_log_columns_match_central_differences_of_the_reference() -> None:
    step = 1e-6
    for time_h in (0.1, 0.5, 0.51, 2.0, 8.0):
        fitted = _concentration_and_log_columns(DOSE_MG, DURATION_H, TRUE_CLEARANCE_L_H, TRUE_VOLUME_L, time_h)
        assert fitted is not None
        concentration, column_cl, column_v = fitted
        assert concentration == pytest.approx(_concentration(time_h), rel=1e-12, abs=1e-12)

        def reference(clearance: float, volume: float, *, time_h: float = time_h) -> float:
            return iv_infusion_concentration(
                dose=DOSE_MG,
                dose_unit="mg",
                infusion_duration=DURATION_H,
                infusion_duration_unit="h",
                volume=volume,
                volume_unit="L",
                clearance=clearance,
                clearance_unit="L/h",
                time=time_h,
                time_unit="h",
                concentration_unit="mg/L",
            )

        d_log_cl = (
            reference(TRUE_CLEARANCE_L_H * math.exp(step), TRUE_VOLUME_L)
            - reference(TRUE_CLEARANCE_L_H * math.exp(-step), TRUE_VOLUME_L)
        ) / (2.0 * step)
        d_log_v = (
            reference(TRUE_CLEARANCE_L_H, TRUE_VOLUME_L * math.exp(step))
            - reference(TRUE_CLEARANCE_L_H, TRUE_VOLUME_L * math.exp(-step))
        ) / (2.0 * step)
        assert column_cl == pytest.approx(d_log_cl, rel=1e-6, abs=1e-8)
        assert column_v == pytest.approx(d_log_v, rel=1e-6, abs=1e-8)


def test_checked_in_infusion_study_omits_blq_and_recovers_clearance() -> None:
    study = PharmacologyStudy.from_files(INFUSION_EXAMPLE_DIR, mode="strict")
    fits = fit_one_compartment_iv_infusion(study, weight="1/y^2")
    assert len(fits) == 2
    expected = {
        subject_id: (body_kg * CLEARANCE_L_PER_H_PER_KG, body_kg * VOLUME_L_PER_KG)
        for subject_id, _sex, body_kg in SUBJECTS
    }
    for fit in fits:
        clearance, volume = expected[fit.subject_id]
        values = _values(fit)
        assert values["CL"] == pytest.approx(clearance, rel=1e-4)
        assert values["V"] == pytest.approx(volume, rel=1e-4)
        assert all(estimate.semantic is not ParameterSemantic.CL_OVER_F for estimate in fit.estimates)
        assert all(
            "T1440" not in observation_id and "PRE" not in observation_id for observation_id in fit.observation_ids
        )
        assert min(fit.observed_mg_l) > 0
    blq = [observation for observation in study.observations if observation.status.value == "BLQ"]
    assert blq
    assert all(observation.numeric_value is None for observation in blq)


def test_evidence_uses_absolute_clearance() -> None:
    fit = fit_one_compartment_iv_infusion_profile(_profile([0.25, 0.5, 1.0, 2.0, 4.0]))
    records = evidence_from_pk_fit(fit, run_id="infusion-fit")
    assert {record.measurement_type for record in records} == {"clearance", "volume"}
    assert {record.claim_type for record in records} == {"CL", "V"}
    assert all(record.numeric_value is not None and record.numeric_value > 0 for record in records)


@pytest.mark.parametrize(
    ("times", "concentrations", "duration", "reason"),
    [
        ([0.25, 0.5], None, DURATION_H, "Fewer than three"),
        ([0.25, 0.5, 1.0], [0.2, 0.4, 0.3], None, "missing or not positive"),
        ([0.25, 0.5, 1.0], [0.2, 0.4, 0.3], 1.0, "not shorter than the observation window"),
        ([0.25, 0.5, 1.0], [0.2, 0.4, 0.3], 2.0, "not shorter than the observation window"),
        ([0.25, 1.0, 2.0, 4.0], [1.0, 1.0, 1.0, 1.0], DURATION_H, "enough curvature"),
        ([0.25, 0.5, 1.0, 2.0, 4.0], [0.4, 0.8, 1.0, 1.2, 1.5], DURATION_H, "not declining"),
    ],
)
def test_failed_curves_do_not_return_clearance(
    times: list[float],
    concentrations: list[float] | None,
    duration: float | None,
    reason: str,
) -> None:
    fit = fit_one_compartment_iv_infusion_profile(
        _profile(times, concentrations, duration_h=duration),
        weight="1/y^2",
    )
    assert fit.estimates == ()
    assert fit.predicted_mg_l == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL
    assert fit.fit.converged is False
    assert any(reason in item for item in fit.fit.limitations)
    records = evidence_from_pk_fit(fit, run_id="infusion-fail")
    assert len(records) == 1
    assert records[0].measurement_type == "clearance"
    assert records[0].numeric_value is None


def test_repeated_regimen_does_not_return_clearance() -> None:
    fit = fit_one_compartment_iv_infusion_profile(_profile([0.25, 0.5, 1.0, 2.0], regimen="repeated_equal"))
    assert fit.estimates == ()
    assert fit.fit.diagnostic_status is DiagnosticStatus.FAIL


def test_bolus_route_is_rejected() -> None:
    with pytest.raises(AnalysisError, match="iv_infusion"):
        fit_one_compartment_iv_infusion_profile(_profile([0.25, 0.5, 1.0, 2.0], route=Route.IV_BOLUS))


def test_unknown_weight_is_rejected() -> None:
    with pytest.raises(AnalysisError, match="weight"):
        fit_one_compartment_iv_infusion_profile(_profile([0.25, 0.5, 1.0, 2.0]), weight="1/sqrt(y)")


def test_cli_fits_the_infusion_example(capsys: pytest.CaptureFixture[str]) -> None:
    assert (
        main(
            [
                "fit",
                str(INFUSION_EXAMPLE_DIR),
                "--model",
                "one_compartment_iv_infusion",
                "--method",
                "weighted_nls",
                "--weight",
                "1/y^2",
            ]
        )
        == 0
    )
    captured = capsys.readouterr().out
    assert "model=one_compartment_iv_infusion" in captured
    assert "CL=" in captured and "V=" in captured
    assert "CL/F" not in captured
    assert (
        main(
            [
                "fit",
                str(INFUSION_EXAMPLE_DIR),
                "--model",
                "one_compartment_iv_infusion",
                "--method",
                "log_linear_ols",
            ]
        )
        == 2
    )
