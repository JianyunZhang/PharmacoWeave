"""Unit tests for multi-scale EphysWeave and OmicsWeave endpoint alignment."""

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.evidence import (
    EvidenceExchange,
    EvidenceQc,
    EvidenceRecord,
    EvidenceType,
    IdentityLink,
    ProjectName,
    ReviewStatus,
)
from pharmacoweave.errors import AlignmentError, ContractViolation
from pharmacoweave.integrations.alignment import (
    align_cross_project_endpoint,
    align_study_cross_project_evidence,
)


def _setup_study_with_identity():
    study = build_study()
    # Add an accepted identity link between EPHYSWEAVE/OMICSWEAVE and SUBJ-01
    link_ephys = IdentityLink(
        link_id="LINK_EW_01",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="STUDY_EPHYS_99",
        source_subject_id="EW_SUBJ_A",
        canonical_subject_id="SUBJ-01",
        match_method="rfid_implant_match",
        review_status=ReviewStatus.ACCEPTED,
    )
    link_omics = IdentityLink(
        link_id="LINK_OW_01",
        source_project=ProjectName.OMICSWEAVE,
        source_study_id="STUDY_OMICS_88",
        source_subject_id="OW_SUBJ_A",
        canonical_subject_id="SUBJ-01",
        match_method="cage_ear_tag_match",
        review_status=ReviewStatus.ACCEPTED,
    )
    link_unreviewed = IdentityLink(
        link_id="LINK_EW_UNREV",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="STUDY_EPHYS_99",
        source_subject_id="EW_SUBJ_UNREV",
        canonical_subject_id="SUBJ-02",
        match_method="manual_guess",
        review_status=ReviewStatus.UNREVIEWED,
    )
    return study.model_copy(
        update={"identity_links": [*study.identity_links, link_ephys, link_omics, link_unreviewed]}
    )


def test_multiscale_alignment_success() -> None:
    study = _setup_study_with_identity()
    ephys_record = EvidenceRecord(
        evidence_id="EVID_EW_001",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ew_42",
        entity_type="lfp_channel",
        entity_id="CH_CA1_THETA",
        measurement_type="lfp_theta_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys_endpoint",
        subject_id="EW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=12.5,
        numeric_value=450.2,
        unit="uV",
        source_artifact_id="artifact://ephys/recordings_session1.nwb",
        limitations=["Research use only."],
    )

    aligned = align_cross_project_endpoint(
        study=study,
        foreign_record=ephys_record,
        dose_event_id="DOSE-SUBJ-01",
        exposure_ref="nca:SUBJ-01:cmax",
        max_allowed_sync_error_ms=100.0,
    )

    assert aligned.canonical_subject_id == "SUBJ-01"
    assert aligned.source_project == ProjectName.EPHYSWEAVE
    assert aligned.source_entity_id == "CH_CA1_THETA"
    assert aligned.source_artifact_id == "artifact://ephys/recordings_session1.nwb"
    assert aligned.sync_error_ms == 12.5
    assert aligned.evidence_grade == "observed"
    assert (
        aligned.traceability_path
        == "dose:DOSE-SUBJ-01 -> exposure:nca:SUBJ-01:cmax -> ephysweave:CH_CA1_THETA -> artifact:artifact://ephys/recordings_session1.nwb"
    )
    assert "dose:DOSE-SUBJ-01" in aligned.derived_from
    assert "exposure:nca:SUBJ-01:cmax" in aligned.derived_from
    assert "artifact:artifact://ephys/recordings_session1.nwb" in aligned.derived_from


def test_multiscale_alignment_rejected_sync_error() -> None:
    study = _setup_study_with_identity()
    ephys_record = EvidenceRecord(
        evidence_id="EVID_EW_SYNC_FAIL",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ew_42",
        entity_type="spike_train",
        entity_id="UNIT_01",
        measurement_type="firing_rate",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys_endpoint",
        subject_id="EW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=850.0,
        numeric_value=15.4,
        unit="1/s",
        source_artifact_id="artifact://ephys/spikes.nwb",
    )

    with pytest.raises(AlignmentError, match="exceeds maximum allowed threshold"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=ephys_record,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
            max_allowed_sync_error_ms=500.0,  # lower than 850.0
        )


def test_multiscale_alignment_rejected_unreviewed_identity() -> None:
    study = _setup_study_with_identity()
    ephys_record = EvidenceRecord(
        evidence_id="EVID_EW_UNREV",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ew_42",
        entity_type="lfp",
        entity_id="CH_01",
        measurement_type="lfp_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys_endpoint",
        subject_id="EW_SUBJ_UNREV",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=5.0,
        numeric_value=12.0,
        unit="uV",
        source_artifact_id="artifact://ephys/ch1.nwb",
    )

    with pytest.raises(AlignmentError, match="requires review_status='accepted'"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=ephys_record,
            dose_event_id="DOSE-SUBJ-02",
            exposure_ref="nca:SUBJ-02:cmax",
        )


def test_multiscale_alignment_rejected_unknown_timebase() -> None:
    study = _setup_study_with_identity()
    ephys_record = EvidenceRecord(
        evidence_id="EVID_EW_BAD_CLOCK",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ew_42",
        entity_type="lfp",
        entity_id="CH_01",
        measurement_type="lfp_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys_endpoint",
        subject_id="EW_SUBJ_A",
        timebase_id="UNKNOWN_UNREGISTERED_CLOCK",
        sync_error_ms=5.0,
        numeric_value=12.0,
        unit="uV",
        source_artifact_id="artifact://ephys/ch1.nwb",
    )

    with pytest.raises(AlignmentError, match="not registered in study timebases"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=ephys_record,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
        )


def test_flux_is_not_concentration_contract_violation() -> None:
    study = _setup_study_with_identity()
    # Omics metabolic flux with concentration units or concentration claim must be rejected!
    flux_with_conc_unit = EvidenceRecord(
        evidence_id="EVID_OW_FLUX_BAD_UNIT",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="run_ow_10",
        entity_type="metabolic_reaction",
        entity_id="RXN_GLYCOLYSIS",
        measurement_type="metabolite_flux",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="flux_rate",
        subject_id="OW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=20.0,
        numeric_value=3.4,
        unit="ug/mL",  # illegal concentration unit for flux!
        source_artifact_id="artifact://omics/flux_model.xml",
    )

    with pytest.raises(ContractViolation, match="Flux is not concentration"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=flux_with_conc_unit,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
        )

    flux_with_conc_claim = EvidenceRecord(
        evidence_id="EVID_OW_FLUX_BAD_CLAIM",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="run_ow_10",
        entity_type="concentration",  # illegal entity_type for flux!
        entity_id="RXN_TCA",
        measurement_type="metabolic_flux",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="concentration",  # illegal claim_type for flux!
        subject_id="OW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=20.0,
        numeric_value=1.5,
        unit="1/h",
        source_artifact_id="artifact://omics/flux_model.xml",
    )

    with pytest.raises(ContractViolation, match="Flux is not concentration"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=flux_with_conc_claim,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
        )


def test_batch_alignment_report() -> None:
    study = _setup_study_with_identity()
    good_record = EvidenceRecord(
        evidence_id="EVID_GOOD",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="run_ow_10",
        entity_type="gene_transcript",
        entity_id="Cyp3a4",
        measurement_type="log2_fold_change",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="differential_expression",
        subject_id="OW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=15.0,
        numeric_value=2.45,
        unit="fold",
        source_artifact_id="artifact://omics/deseq2_results.tsv",
    )
    bad_sync_record = EvidenceRecord(
        evidence_id="EVID_BAD_SYNC",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ew_10",
        entity_type="lfp",
        entity_id="CH_01",
        measurement_type="lfp_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys_endpoint",
        subject_id="EW_SUBJ_A",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=5000.0,  # exceeds max 1000.0
        numeric_value=100.0,
        unit="uV",
        source_artifact_id="artifact://ephys/rec.nwb",
    )

    exchange = EvidenceExchange(records=[good_record, bad_sync_record])
    mapping = {
        "EVID_GOOD": ("DOSE-SUBJ-01", "nca:SUBJ-01:auc_inf"),
        "EVID_BAD_SYNC": ("DOSE-SUBJ-01", "nca:SUBJ-01:cmax"),
    }

    report = align_study_cross_project_evidence(
        study=study,
        exchange=exchange,
        dose_exposure_mapping=mapping,
        max_allowed_sync_error_ms=1000.0,
    )

    assert len(report.aligned_endpoints) == 1
    assert report.aligned_endpoints[0].source_entity_id == "Cyp3a4"
    assert len(report.rejected_alignments) == 1
    assert "exceeds maximum allowed threshold" in report.rejected_alignments[0].reason
