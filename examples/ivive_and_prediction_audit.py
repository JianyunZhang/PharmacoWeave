"""Exploratory Well-Stirred IVIVE and Prior-vs-Estimated Prediction Audit.

Phase 2 exploratory capability demonstration:
1. In vitro ADME records (CLint, fu_plasma, fu_mic, blood-to-plasma ratio).
2. Explicit-assumption exploratory well-stirred hepatic IVIVE (stays predicted).
3. Calibration audit comparing priors (TheraWeave CADD + IVIVE) against in vivo NCA estimates.
4. Packaging into an immutable ResultBundle with HTML report and diagnostic tables.
"""

from pathlib import Path

from examples.iv_bolus_case import build_study
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.diagnostics.prediction_audit import audit_study_predictions
from pharmacoweave.dmpk import (
    AdmeAssayKind,
    AdmeStatus,
    InVitroAdmeRecord,
    PhysiologicalScalingParams,
    predict_hepatic_clearance_well_stirred,
)
from pharmacoweave.integrations.guards import load_exchange
from pharmacoweave.integrations.pharmacology_evidence import collect_run_evidence, evidence_from_ivive

REPO_ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    # 1. Load synthetic IV bolus study
    study = build_study()
    print(f"Loaded study: {study.metadata.study_id} ({study.metadata.species})")

    # 2. In vitro ADME records
    clint_rec = InVitroAdmeRecord(
        record_id="ADME-MIC-RAT",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="rat liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=[RESEARCH_USE_LIMITATION],
        value=28.5,
        unit="uL/min/mg_protein",
    )
    fu_p_rec = InVitroAdmeRecord(
        record_id="ADME-FUP-RAT",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.FU_PLASMA,
        system="rat plasma",
        status=AdmeStatus.MEASURED,
        limitations=[RESEARCH_USE_LIMITATION],
        value=12.0,
        unit="percent",
    )
    fu_mic_rec = InVitroAdmeRecord(
        record_id="ADME-FUMIC-RAT",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.FU_MIC,
        system="rat liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=[RESEARCH_USE_LIMITATION],
        value=0.85,
        unit="fraction",
    )
    bp_rec = InVitroAdmeRecord(
        record_id="ADME-BP-RAT",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO,
        system="rat whole blood",
        status=AdmeStatus.MEASURED,
        limitations=[RESEARCH_USE_LIMITATION],
        value=1.1,
        unit="ratio",
    )

    # 3. Caller-supplied physiology for rat (Wistar)
    rat_phys = PhysiologicalScalingParams(
        profile_id="PHYS-RAT-WISTAR",
        species="Rattus norvegicus",
        body_weight_kg=0.25,
        liver_weight_g=10.0,
        hepatic_blood_flow_l_h=1.05,
        microsomal_protein_mg_per_g_liver=45.0,
        source="Davies & Morris (1993), Pharm Res 10(7):1093-1095",
    )

    # 4. Exploratory IVIVE
    ivive_pred = predict_hepatic_clearance_well_stirred(
        clint_rec,
        physiology=rat_phys,
        fu_plasma=fu_p_rec,
        fu_inc=fu_mic_rec,
        blood_to_plasma_ratio=bp_rec,
    )
    print("\n--- Exploratory IVIVE Prediction (Well-Stirred) ---")
    print(f"Prediction ID:             {ivive_pred.prediction_id}")
    print(f"Evidence Grade:            {ivive_pred.evidence_grade} (never in vivo estimated/observed)")
    print(f"CLint,u,vitro:             {ivive_pred.clint_u_vitro_value:.4f} uL/min/mg_protein")
    print(f"CLint,u,H:                 {ivive_pred.clint_u_hepatic_l_h:.4f} L/h")
    print(f"Hepatic Extraction Ratio:  {ivive_pred.hepatic_extraction_ratio:.4f}")
    print(f"Hepatic Blood Clearance:   {ivive_pred.hepatic_blood_clearance_l_h:.4f} L/h")
    print(f"Hepatic Plasma Clearance:  {ivive_pred.hepatic_plasma_clearance_l_h:.4f} L/h")
    print(f"Hepatic Plasma CL/kg:      {ivive_pred.hepatic_plasma_clearance_l_h_kg:.4f} L/h/kg")

    # 5. In vivo single-dose NCA
    nca_config_path = REPO_ROOT / "configs" / "nca_single_dose_iv.yaml"
    nca_results = study.run_nca(nca_config_path)
    in_vivo_estimates = collect_run_evidence(run_id="run-nca", nca_parameters=nca_results)

    # 6. Prior-vs-Estimated Prediction Audit
    # Load quantitative priors from TheraWeave CADD demo
    priors_exchange = load_exchange(
        REPO_ROOT / "examples" / "cross_project_demo" / "theraweave_quantitative_priors.json"
    )
    priors = list(priors_exchange.records)
    # Also attach the exploratory IVIVE prediction as a prior
    priors.append(evidence_from_ivive(ivive_pred, run_id="run-ivive"))

    audit_report = audit_study_predictions(
        study,
        in_vivo_estimates,
        priors=priors,
        audit_id="AUDIT-DEMO-001",
    )

    print("\n--- Prediction Audit Summary Diagnostics ---")
    for s in audit_report.summaries:
        print(
            f"Measurement: {s.measurement_type:<15} N={s.n_comparisons} "
            f"AFE={s.afe:.3f} AAFE={s.aafe:.3f} "
            f"Within 2-fold: {s.fraction_within_2_fold * 100:.1f}% "
            f"Within 3-fold: {s.fraction_within_3_fold * 100:.1f}% "
            f"Outliers: {s.n_outlier_2_fold}"
        )

    print("\n--- Audited Comparisons ---")
    for e in audit_report.entries:
        print(
            f"{e.prior_evidence_id} vs {e.estimate_evidence_id} ({e.subject_id}): "
            f"Predicted={e.predicted_value:.4f} {e.unit}, Estimated={e.estimated_value:.4f} {e.unit}, "
            f"Ratio={e.ratio:.3f}, FoldError={e.fold_error:.3f} [{e.bias_direction.value}]"
        )


if __name__ == "__main__":
    main()
