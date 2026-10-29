"""Fits that must not be marked pass."""

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.pk_models import fit_one_compartment_iv_bolus


def test_one_sparse_profile_is_not_identifiable_and_is_not_a_pass() -> None:
    study = build_study()
    samples = [
        sample
        for sample in study.biosamples
        if sample.subject_id != "SUBJ-02" or sample.time_since_dose in {-5, 5}
    ]
    identifiers = {sample.sample_id for sample in samples}
    observations = [item for item in study.observations if item.sample_id in identifiers]
    thinned = study.model_copy(update={"biosamples": samples, "observations": observations})
    fitted = fit_one_compartment_iv_bolus(thinned)
    failed = [item for item in fitted.records.model_fits if "SUBJ-02" in item.fit_id]
    passed = [item for item in fitted.records.model_fits if "SUBJ-01" in item.fit_id]
    assert failed[0].diagnostic_status is DiagnosticStatus.FAIL
    assert failed[0].diagnostic_status is not DiagnosticStatus.PASS
    assert any("not identifiable" in item for item in failed[0].limitations)
    assert not any(item.fit_id == failed[0].fit_id for item in fitted.records.parameter_estimates)
    assert not any(item.fit_id == failed[0].fit_id for item in fitted.predictions)
    assert passed[0].diagnostic_status is DiagnosticStatus.PASS
    assert any(item.fit_id == passed[0].fit_id for item in fitted.records.parameter_estimates)


def test_optimizer_that_cannot_finish_is_not_a_pass() -> None:
    fitted = fit_one_compartment_iv_bolus(build_study(), initial=(1e-3, 1e-3), max_nfev=1)
    assert fitted.records.model_fits
    assert all(item.diagnostic_status is DiagnosticStatus.FAIL for item in fitted.records.model_fits)
    assert all(item.converged is False for item in fitted.records.model_fits)
    assert fitted.records.parameter_estimates == []
    assert fitted.predictions == ()


def test_flat_profile_does_not_pass() -> None:
    study = build_study()
    observations = []
    for observation in study.observations:
        if observation.status is ObservationStatus.QUANTIFIED:
            observations.append(
                observation.model_copy(update={"numeric_value": 1.0, "raw_value": "1.0", "reported_value": "1.0"})
            )
        else:
            observations.append(observation)
    fitted = fit_one_compartment_iv_bolus(study.model_copy(update={"observations": observations}))
    assert fitted.records.model_fits
    assert all(item.diagnostic_status is not DiagnosticStatus.PASS for item in fitted.records.model_fits)
    assert fitted.records.parameter_estimates == []


def test_initial_must_be_positive() -> None:
    with pytest.raises(ValueError, match="positive"):
        fit_one_compartment_iv_bolus(build_study(), initial=(0.0, 1.0))
