"""Multi-scale electrophysiology and omics alignment demonstration.

Demonstrates:
1. Loading an in vivo PK study with reviewed identity mappings.
2. Computing NCA exposures (Cmax, AUC).
3. Aligning EphysWeave LFP endpoints and OmicsWeave transcriptomic responses.
4. Verifying clock synchrony tolerance (sync_error_ms) and provenance lineage.
5. Verifying negative boundary: refusing to classify Omics metabolic flux as concentration.
"""

from pathlib import Path

from examples.iv_bolus_case import build_study
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


def run_demo() -> None:
    print("=== PharmacoWeave 0.1.0: Multi-Scale Time Alignment Demo ===")

    # 1. Prepare in vivo study with identity map links
    study = build_study()
    ephys_link = IdentityLink(
        link_id="LINK_EW_001",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="EPHYS_RIG_A",
        source_subject_id="EW_ANIMAL_01",
        canonical_subject_id="SUBJ-01",
        match_method="rfid_chip_match",
        review_status=ReviewStatus.ACCEPTED,
        notes="Verified identity between PK catheterization and neural telemetry.",
    )
    omics_link = IdentityLink(
        link_id="LINK_OW_001",
        source_project=ProjectName.OMICSWEAVE,
        source_study_id="OMICS_FLOWCELL_B",
        source_subject_id="OW_ANIMAL_01",
        canonical_subject_id="SUBJ-01",
        match_method="ear_biopsy_and_cage",
        review_status=ReviewStatus.ACCEPTED,
        notes="Verified liver sample taken post-necropsy from SUBJ-01.",
    )
    study = study.model_copy(
        update={"identity_links": [*study.identity_links, ephys_link, omics_link]}
    )
    print(f"Loaded study: {study.metadata.study_id} with {len(study.identity_links)} identity links.")

    # 2. Compute in vivo exposure via NCA
    config_path = Path(__file__).resolve().parents[1] / "configs" / "nca_single_dose_iv.yaml"
    nca_results = study.run_nca(config_path)
    subj_01_nca = next(r for r in nca_results if r.subject_id == "SUBJ-01")
    print(f"Subject SUBJ-01 Exposure: Cmax={subj_01_nca.cmax:.4f} {subj_01_nca.cmax_unit}, "
          f"AUCinf={subj_01_nca.auc_inf:.4f} {subj_01_nca.auc_inf_unit}")

    # 3. Create foreign evidence from EphysWeave and OmicsWeave
    ephys_endpoint = EvidenceRecord(
        evidence_id="EVID_EPHYS_THETA_01",
        source_project=ProjectName.EPHYSWEAVE,
        source_run_id="session_rec_20261031_01",
        entity_type="lfp_band",
        entity_id="CA1_THETA",
        measurement_type="lfp_theta_power",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="electrophysiology_endpoint",
        subject_id="EW_ANIMAL_01",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=14.2,  # 14.2 ms synchronization jitter
        numeric_value=284.5,
        unit="uV",
        source_artifact_id="sha256:d41d8cd98f00b204e9800998ecf8427e",
        limitations=["Research use only."],
    )

    omics_endpoint = EvidenceRecord(
        evidence_id="EVID_OMICS_CYP_01",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="deseq2_run_20261031_02",
        entity_type="mrna_transcript",
        entity_id="CYP3A4",
        measurement_type="gene_expression_log2fc",
        evidence_type=EvidenceType.OBSERVED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="differential_expression",
        subject_id="OW_ANIMAL_01",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=45.0,  # 45 ms clock uncertainty
        numeric_value=1.92,
        unit="fold",
        source_artifact_id="sha256:a69f73ccd356db7c5e819d7fa071b281",
        limitations=["Research use only."],
    )

    # 4. Perform cross-scale alignment
    aligned_ephys = align_cross_project_endpoint(
        study=study,
        foreign_record=ephys_endpoint,
        dose_event_id="DOSE-SUBJ-01",
        exposure_ref="nca:SUBJ-01:cmax",
        max_allowed_sync_error_ms=50.0,
    )
    print("\nAligned Ephys Endpoint:")
    print(f"  Aligned ID:        {aligned_ephys.aligned_id}")
    print(f"  Canonical Subject: {aligned_ephys.canonical_subject_id}")
    print(f"  Traceability:      {aligned_ephys.traceability_path}")
    print(f"  Sync Error:        {aligned_ephys.sync_error_ms} ms")

    aligned_omics = align_cross_project_endpoint(
        study=study,
        foreign_record=omics_endpoint,
        dose_event_id="DOSE-SUBJ-01",
        exposure_ref="nca:SUBJ-01:auc_inf",
        max_allowed_sync_error_ms=50.0,
    )
    print("\nAligned Omics Endpoint:")
    print(f"  Aligned ID:        {aligned_omics.aligned_id}")
    print(f"  Canonical Subject: {aligned_omics.canonical_subject_id}")
    print(f"  Traceability:      {aligned_omics.traceability_path}")
    print(f"  Sync Error:        {aligned_omics.sync_error_ms} ms")

    # 5. Verify clock sync limit rejection
    jittery_ephys = ephys_endpoint.model_copy(update={"sync_error_ms": 1200.0})
    try:
        align_cross_project_endpoint(
            study=study,
            foreign_record=jittery_ephys,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
            max_allowed_sync_error_ms=500.0,
        )
        raise RuntimeError("Expected AlignmentError on excessive jitter!")
    except AlignmentError as exc:
        print(f"\n[Guardrail Verified] Clock sync limit exceeded: {exc}")

    # 6. Verify negative constraint: Omics metabolic flux is NOT concentration
    illegal_flux = EvidenceRecord(
        evidence_id="EVID_OW_FLUX_01",
        source_project=ProjectName.OMICSWEAVE,
        source_run_id="fba_model_01",
        entity_type="concentration",  # ILLEGAL
        entity_id="GLUCOSE_UPTAKE",
        measurement_type="metabolite_flux",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.ACCEPTED,
        claim_type="concentration",  # ILLEGAL
        subject_id="OW_ANIMAL_01",
        timebase_id="STUDY-CLOCK",
        sync_error_ms=10.0,
        numeric_value=2.5,
        unit="umol/L",  # ILLEGAL FOR FLUX
        source_artifact_id="sha256:11223344",
    )
    try:
        align_cross_project_endpoint(
            study=study,
            foreign_record=illegal_flux,
            dose_event_id="DOSE-SUBJ-01",
            exposure_ref="nca:SUBJ-01:cmax",
        )
        raise RuntimeError("Expected ContractViolation on flux marked as concentration!")
    except ContractViolation as exc:
        print(f"[Guardrail Verified] Negative constraint: {exc}")

    print("\nAll multi-scale alignment validations passed successfully.")


if __name__ == "__main__":
    run_demo()
