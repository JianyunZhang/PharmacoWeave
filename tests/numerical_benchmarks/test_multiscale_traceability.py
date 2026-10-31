"""Analytical benchmarks for Phase 3 multi-scale time alignment and end-to-end traceability.

Verifies mathematical and boundary contracts for:
1. Dose event -> exposure -> EphysWeave/OmicsWeave PD endpoint -> source artifact.
2. Explicit clock timebase and synchronization error (sync_error_ms) thresholding.
3. Reviewed identity link enforcement (review_status='accepted').
4. Negative boundary: Omics metabolic flux is not concentration.
5. Study observations bridge to EffectObservation and ToxEndpoint.
"""

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import (
    Observation,
    ObservationStatus,
    ObservationType,
)
from pharmacoweave.core.evidence import (
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
)
from pharmacoweave.pd_models.records import EffectDimension, EffectTimeScale
from pharmacoweave.toxicology.records import ToxKind, ToxSeverity


def test_analytical_multiscale_traceability_chain() -> None:
    study = build_study()
    # Add identity links with reviewed accepted status
    link_ephys = IdentityLink(
        link_id="ID-LINK-EPHYS-001",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="STUDY-EPHYS-CA1",
        source_subject_id="EW-RAT-42",
        canonical_subject_id="SUBJ-01",
        match_method="rfid_and_tattoo",
        review_status=ReviewStatus.ACCEPTED,
        notes="Verified identical animal across PK and electrophysiology rigs.",
    )
    link_omics = IdentityLink(
        link_id="ID-LINK-OMICS-001",
        source_project=ProjectName.OMICSWEAVE,
        source_study_id="STUDY-OMICS-LIVER",
        source_subject_id="OW-RAT-42",
        canonical_subject_id="SUBJ-01",
        match_method="ear_notch_and_cage",
        review_status=ReviewStatus.ACCEPTED,
        notes="Verified liver biopsy from same subject.",
    )
    study = study.model_copy(
        update={"identity_links": [*study.identity_links, link_ephys, link_omics]}
    )

    # 1. EphysWeave neural endpoint
    # Ephys measures hippocampal CA1 theta power after administration
    ephys_record = EvidenceRecord(
        evidence_id="EW-EVID-THETA-01",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_ephys_theta_session_01",
        entity_type="lfp_power_spectrum",
        entity_id="CA1_THETA_BAND",
        measurement_type="lfp_theta_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="electrophysiology_endpoint",
        subject_id="EW-RAT-42",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=8.45,  # within 50.0 ms tolerance
        numeric_value=314.159,
        unit="uV",
        source_artifact_id="sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        limitations=["Research use only."],
    )

    # 2. OmicsWeave transcriptomic response
    # Omics measures Cyp3a4 induction log2 fold change
    omics_record = EvidenceRecord(
        evidence_id="OW-EVID-CYP-01",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="run_rnaseq_deseq2_01",
        entity_type="mrna_transcript",
        entity_id="CYP3A4",
        measurement_type="gene_expression_log2fc",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="differential_gene_expression",
        subject_id="OW-RAT-42",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=32.10,  # within 50.0 ms tolerance
        numeric_value=1.85,
        unit="fold",
        source_artifact_id="sha256:4b227777d4dd1fc61c6f884f48641d02b4d121d3fd328cb08b5531fcacdabf8a",
        limitations=["Research use only."],
    )

    # Align Ephys endpoint to dose event DOSE-SUBJ-01 and exposure Cmax
    aligned_ephys = align_cross_project_endpoint(
        study=study,
        foreign_record=ephys_record,
        dose_event_id="DOSE-SUBJ-01",
        exposure_ref="nca:SUBJ-01:cmax",
        max_allowed_sync_error_ms=50.0,
    )

    # Verify exact lineage and properties
    assert aligned_ephys.canonical_subject_id == "SUBJ-01"
    assert aligned_ephys.dose_event_id == "DOSE-SUBJ-01"
    assert aligned_ephys.exposure_ref == "nca:SUBJ-01:cmax"
    assert aligned_ephys.sync_error_ms == 8.45
    assert aligned_ephys.numeric_value == pytest.approx(314.159, rel=1e-6)
    assert aligned_ephys.unit == "uV"
    assert aligned_ephys.evidence_grade == "observed"

    expected_ephys_path = (
        "dose:DOSE-SUBJ-01 -> exposure:nca:SUBJ-01:cmax -> "
        "ephysweave:CA1_THETA_BAND -> "
        "artifact:sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    )
    assert aligned_ephys.traceability_path == expected_ephys_path
    assert aligned_ephys.derived_from == [
        "dose:DOSE-SUBJ-01",
        "exposure:nca:SUBJ-01:cmax",
        "ephysweave:CA1_THETA_BAND",
        "artifact:sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    ]

    # Align Omics endpoint to dose event DOSE-SUBJ-01 and exposure AUC_inf
    aligned_omics = align_cross_project_endpoint(
        study=study,
        foreign_record=omics_record,
        dose_event_id="DOSE-SUBJ-01",
        exposure_ref="nca:SUBJ-01:auc_inf",
        max_allowed_sync_error_ms=50.0,
    )
    assert aligned_omics.canonical_subject_id == "SUBJ-01"
    assert aligned_omics.sync_error_ms == 32.10
    assert aligned_omics.numeric_value == pytest.approx(1.85, rel=1e-6)

    expected_omics_path = (
        "dose:DOSE-SUBJ-01 -> exposure:nca:SUBJ-01:auc_inf -> "
        "omicsweave:CYP3A4 -> "
        "artifact:sha256:4b227777d4dd1fc61c6f884f48641d02b4d121d3fd328cb08b5531fcacdabf8a"
    )
    assert aligned_omics.traceability_path == expected_omics_path


def test_analytical_guardrails_sync_and_identity() -> None:
    study = build_study()
    # 1. Unreviewed identity link must be rejected
    link_unreviewed = IdentityLink(
        link_id="ID-LINK-UNREVIEWED",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="STUDY-EPHYS-CA1",
        source_subject_id="EW-RAT-UNREVIEWED",
        canonical_subject_id="SUBJ-01",
        match_method="automated_guess",
        review_status=ReviewStatus.UNREVIEWED,
    )
    study = study.model_copy(update={"identity_links": [*study.identity_links, link_unreviewed]})

    ephys_rec = EvidenceRecord(
        evidence_id="EW-REC-01",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="run_01",
        entity_type="lfp",
        entity_id="THETA",
        measurement_type="lfp_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="ephys",
        subject_id="EW-RAT-UNREVIEWED",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=5.0,
        numeric_value=100.0,
        unit="uV",
        source_artifact_id="sha256:deadbeef",
    )

    with pytest.raises(AlignmentError, match="requires review_status='accepted'"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=ephys_rec,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="cmax",
        )

    # 2. Clock synchrony error over limit
    link_accepted = IdentityLink(
        link_id="ID-LINK-ACCEPTED",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="STUDY-EPHYS-CA1",
        source_subject_id="EW-RAT-ACCEPTED",
        canonical_subject_id="SUBJ-01",
        match_method="verified",
        review_status=ReviewStatus.ACCEPTED,
    )
    study = study.model_copy(update={"identity_links": [*study.identity_links, link_accepted]})

    ephys_rec_high_jitter = ephys_rec.model_copy(
        update={"subject_id": "EW-RAT-ACCEPTED", "sync_error_ms": 120.5}
    )
    with pytest.raises(AlignmentError, match="Clock synchronization error 120.50 ms exceeds"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=ephys_rec_high_jitter,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="cmax",
            max_allowed_sync_error_ms=100.0,
        )


def test_analytical_negative_constraint_flux_not_concentration() -> None:
    study = build_study()
    link = IdentityLink(
        link_id="ID-LINK-OW",
        source_project=ProjectName.OMICSWEAVE,
        source_study_id="STUDY-OW",
        source_subject_id="OW-SUBJ-01",
        canonical_subject_id="SUBJ-01",
        match_method="verified",
        review_status=ReviewStatus.ACCEPTED,
    )
    study = study.model_copy(update={"identity_links": [*study.identity_links, link]})

    # Metabolic flux claimed as concentration
    flux_as_conc = EvidenceRecord(
        evidence_id="OW-FLUX-ILLEGAL-CONC",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="run_fba",
        entity_type="concentration",  # ILLEGAL
        entity_id="GLUCOSE_FLUX",
        measurement_type="metabolite_flux",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="concentration",  # ILLEGAL
        subject_id="OW-SUBJ-01",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=10.0,
        numeric_value=50.0,
        unit="umol/L",  # ILLEGAL FOR FLUX
        source_artifact_id="sha256:feedbeef",
    )

    with pytest.raises(ContractViolation, match="Flux is not concentration"):
        align_cross_project_endpoint(
            study=study,
            foreign_record=flux_as_conc,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="auc",
        )


def test_analytical_observation_bridge_exact_preservation() -> None:
    study = build_study()
    # Add an observed PD effect (e.g. 72.5% inhibition) and toxicology endpoint
    eff_obs = Observation(
        observation_id="OBS-EFF-100",
        sample_id="SUBJ-01-T0015",
        observation_type=ObservationType.EFFECT,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=72.5,
        unit="percent",
    )
    tox_obs = Observation(
        observation_id="OBS-TOX-200",
        sample_id="SUBJ-01-T0060",
        observation_type=ObservationType.TOXICITY,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=14.2,
        unit="percent",
        reported_value="Transient body weight loss",
    )
    study = study.model_copy(update={"observations": [*study.observations, eff_obs, tox_obs]})

    # Bridge effect observations
    effects = study.bridge_effects(
        endpoint_overrides={
            "OBS-EFF-100": {
                "endpoint_name": "TargetInhibition",
            }
        }
    )
    assert len(effects) == 1
    e = effects[0]
    assert e.effect_id == "OBS-EFF-100"
    assert e.endpoint_name == "TargetInhibition"
    assert e.numeric_value == pytest.approx(72.5, rel=1e-6)
    assert e.unit == "percent"
    assert e.dimension == EffectDimension.PERCENT
    assert e.time_scale == EffectTimeScale.DIRECT
    assert e.exposure_ref == "DOSE-SUBJ-01"
    assert e.time_since_exposure == 15.0
    assert e.time_since_exposure_unit == "min"

    # Bridge toxicology endpoints
    tox_list = study.bridge_toxicology(
        endpoint_overrides={
            "OBS-TOX-200": {
                "kind": ToxKind.BODY_WEIGHT_CHANGE,
                "severity": ToxSeverity.MILD,
            }
        }
    )
    assert len(tox_list) == 1
    t = tox_list[0]
    assert t.endpoint_id == "OBS-TOX-200"
    assert t.kind == ToxKind.BODY_WEIGHT_CHANGE
    assert t.severity == ToxSeverity.MILD
    assert t.numeric_value == pytest.approx(14.2, rel=1e-6)
    assert t.unit == "percent"
    assert t.compound_id == "PW-0001"
    assert t.dose_event_id == "DOSE-SUBJ-01"
    assert t.subject_id == "SUBJ-01"
