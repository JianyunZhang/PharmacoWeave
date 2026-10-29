"""CL and V recovery for the one-compartment IV bolus fit."""

import math
import random

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.pk_models import fit_one_compartment_iv_bolus, iv_bolus_concentration

# Noise-free concentrations are rounded to 6 decimals. This is the acceptance
# window versus the generating CL and V, not a claim of extra precision.
_NOISE_FREE_REL = 1e-4
# Proportional noise at sigma 0.02 with seed 20261029. The window is wide
# enough for that one draw and is not a population confidence interval.
_NOISY_REL = 0.10


def test_noise_free_study_recovers_clearance_and_volume() -> None:
    fitted = fit_one_compartment_iv_bolus(build_study())
    assert all(item.diagnostic_status is DiagnosticStatus.PASS for item in fitted.records.model_fits)
    assert all(item.converged is True for item in fitted.records.model_fits)
    assert all(item.limitations for item in fitted.records.diagnostics)
    expected = {"SUBJ-01": (0.05, 0.25), "SUBJ-02": (0.06, 0.30)}
    found: dict[str, dict[str, float]] = {}
    for estimate in fitted.records.parameter_estimates:
        subject = "SUBJ-01" if "SUBJ-01" in estimate.fit_id else "SUBJ-02"
        found.setdefault(subject, {})[estimate.parameter] = estimate.value
        assert estimate.semantic.value in {"CL", "V"}
        assert estimate.bound_hit is False
    for subject, (clearance, volume) in expected.items():
        assert found[subject]["CL"] == pytest.approx(clearance, rel=_NOISE_FREE_REL)
        assert found[subject]["V"] == pytest.approx(volume, rel=_NOISE_FREE_REL)
    assert fitted.predictions
    sample = fitted.predictions[0]
    assert sample.predicted == pytest.approx(
        iv_bolus_concentration(
            dose=0.25 if sample.subject_id == "SUBJ-01" else 0.3,
            dose_unit="mg",
            volume=found[sample.subject_id]["V"],
            volume_unit="L",
            clearance=found[sample.subject_id]["CL"],
            clearance_unit="L/h",
            time=sample.time_h,
            time_unit="h",
            concentration_unit=sample.unit,
        ),
        rel=1e-8,
    )


def test_proportional_noise_stays_near_the_generating_parameters() -> None:
    study = build_study()
    generator = random.Random(20261029)
    observations = []
    for observation in study.observations:
        if observation.status is not ObservationStatus.QUANTIFIED or observation.numeric_value is None:
            observations.append(observation)
            continue
        noisy = observation.numeric_value * math.exp(generator.gauss(0.0, 0.02))
        observations.append(
            observation.model_copy(
                update={
                    "numeric_value": noisy,
                    "raw_value": f"{noisy:.6f}",
                    "reported_value": f"{noisy:.6f}",
                }
            )
        )
    fitted = fit_one_compartment_iv_bolus(study.model_copy(update={"observations": observations}))
    by_subject = {"SUBJ-01": (0.05, 0.25), "SUBJ-02": (0.06, 0.30)}
    for estimate in fitted.records.parameter_estimates:
        subject = "SUBJ-01" if "SUBJ-01" in estimate.fit_id else "SUBJ-02"
        expected = by_subject[subject][0] if estimate.parameter == "CL" else by_subject[subject][1]
        assert estimate.value == pytest.approx(expected, rel=_NOISY_REL)
    assert all(item.diagnostic_status is DiagnosticStatus.PASS for item in fitted.records.model_fits)
