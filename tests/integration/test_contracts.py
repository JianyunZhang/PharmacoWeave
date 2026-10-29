"""Cross-project evidence stays at its original grade."""

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.errors import ContractViolation
from pharmacoweave.integrations import (
    attach_biomarker_evidence,
    attach_candidate_bundle,
    attach_neural_endpoints,
    read_neural_endpoint,
    run_backend,
    supported_backends,
)
from pharmacoweave.integrations.guards import load_exchange, promote_to_observation


def test_predicted_admet_does_not_become_a_clearance_or_concentration(repo_root) -> None:
    study = build_study()
    path = repo_root / "examples" / "cross_project_demo" / "theraweave_predicted_admet.json"
    attached = attach_candidate_bundle(study, path)
    assert len(attached.observations) == len(study.observations)
    assert attached.compounds == study.compounds
    assert len(attached.external_evidence) == 1
    record = attached.external_evidence[0]
    assert record.evidence_type.value == "predicted"
    assert record.measurement_type == "clearance"
    with pytest.raises(ContractViolation):
        promote_to_observation(record)


def test_omics_flux_can_be_cited_and_a_mislabeled_concentration_is_rejected(repo_root) -> None:
    from pharmacoweave.integrations.guards import screen_exchange

    study = build_study()
    demo = repo_root / "examples" / "cross_project_demo"
    attached = attach_biomarker_evidence(study, demo / "omics_flux.json")
    assert attached.external_evidence[0].measurement_type == "flux"
    exchange = load_exchange(demo / "omics_mislabeled_concentration.json")
    with pytest.raises(ContractViolation, match="OmicsWeave concentration"):
        screen_exchange(exchange)


def test_ephys_endpoint_requires_sync_and_does_not_align_clocks(repo_root) -> None:
    study = build_study()
    demo = repo_root / "examples" / "cross_project_demo"
    with pytest.raises(ContractViolation, match="sync_error_ms"):
        attach_neural_endpoints(study, demo / "ephys_missing_sync.json")
    attached = attach_neural_endpoints(study, demo / "ephys_endpoint.json")
    endpoint = read_neural_endpoint(attached.external_evidence[0])
    assert endpoint.sync_error_ms == 2.5
    assert endpoint.timebase_id == "PLEXON-CLOCK"
    assert attached.identity_links == study.identity_links
    assert [item.sample_id for item in attached.biosamples] == [item.sample_id for item in study.biosamples]


def test_no_external_engine_is_registered() -> None:
    from pharmacoweave.errors import PlannedCapabilityError

    assert supported_backends() == []
    with pytest.raises(PlannedCapabilityError):
        run_backend("nlmixr2")
