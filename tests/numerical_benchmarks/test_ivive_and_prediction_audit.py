"""Numerical benchmarks for exploratory well-stirred IVIVE and prior prediction audit."""

import math

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.evidence import EvidenceQc, EvidenceRecord, EvidenceType, ProjectName
from pharmacoweave.diagnostics.prediction_audit import (
    BiasDirection,
    PredictionAuditEntry,
    audit_evidence_predictions,
    audit_study_predictions,
    compute_fold_metrics,
    summarize_audit_entries,
)
from pharmacoweave.dmpk import (
    AdmeAssayKind,
    AdmeStatus,
    InVitroAdmeRecord,
    PhysiologicalScalingParams,
    predict_hepatic_clearance_well_stirred,
)
from pharmacoweave.errors import ContractViolation
from pharmacoweave.integrations.pharmacology_evidence import evidence_from_nca


def test_well_stirred_ivive_hand_calculated_benchmark() -> None:
    """Compare well-stirred IVIVE results against manual analytical formulas."""
    # Input parameters
    clint_app = 25.0  # uL/min/mg_protein
    fu_inc = 0.8
    mppgl = 45.0  # mg/g
    liver_weight_g = 1800.0  # g
    fu_p = 0.05
    bp_ratio = 1.25
    qh = 90.0  # L/h
    bw_kg = 70.0  # kg

    # Analytical steps:
    # 1. CLint,u,vitro = 25.0 / 0.8 = 31.25 uL/min/mg_protein
    # 2. CLint,u,H = 31.25 * 45.0 * 1800.0 * 60 / 1e6 = 151.875 L/h
    # 3. fu_b = fu_p / bp_ratio = 0.05 / 1.25 = 0.04
    # 4. CLH,b = QH * fu_b * CLint,u,H / (QH + fu_b * CLint,u,H)
    #    Numerator = 90.0 * 0.04 * 151.875 = 546.75
    #    Denominator = 90.0 + (0.04 * 151.875) = 90.0 + 6.075 = 96.075
    #    CLH,b = 546.75 / 96.075 = 5.6908665105386415 L/h
    # 5. ER = CLH,b / QH = 5.6908665105386415 / 90.0 = 0.06323185011709602
    # 6. CLH,p = bp_ratio * CLH,b = 1.25 * 5.6908665105386415 = 7.113583138173302 L/h
    # 7. CLH,p,kg = CLH,p / 70.0 = 0.1016226162596186 L/h/kg

    clint_rec = InVitroAdmeRecord(
        record_id="ADME-MIC-BENCH",
        compound_id="PW-TEST",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=["Research use only. Not a dosing recommendation."],
        value=clint_app,
        unit="uL/min/mg_protein",
    )
    phys = PhysiologicalScalingParams(
        profile_id="PHYS-BENCH",
        species="Homo sapiens",
        body_weight_kg=bw_kg,
        liver_weight_g=liver_weight_g,
        hepatic_blood_flow_l_h=qh,
        microsomal_protein_mg_per_g_liver=mppgl,
        source="Analytical benchmark test fixture",
    )

    pred = predict_hepatic_clearance_well_stirred(
        clint_rec,
        physiology=phys,
        fu_plasma=fu_p,
        fu_inc=fu_inc,
        blood_to_plasma_ratio=bp_ratio,
    )

    assert pred.clint_u_vitro_value == pytest.approx(31.25, rel=1e-9)
    assert pred.clint_u_hepatic_l_h == pytest.approx(151.875, rel=1e-9)
    assert pred.fu_blood == pytest.approx(0.04, rel=1e-9)
    assert pred.hepatic_blood_clearance_l_h == pytest.approx(5.6908665, rel=1e-6)
    assert pred.hepatic_extraction_ratio == pytest.approx(0.06323185, rel=1e-6)
    assert pred.hepatic_plasma_clearance_l_h == pytest.approx(7.1135831, rel=1e-6)
    assert pred.hepatic_plasma_clearance_l_h_kg == pytest.approx(0.1016226, rel=1e-6)


def test_well_stirred_ivive_extraction_limits() -> None:
    """Verify low and high extraction ratio asymptotic limits."""
    # Low extraction: CLH,b -> fu_b * CLint,u,H
    phys = PhysiologicalScalingParams(
        profile_id="PHYS-HUMAN",
        species="Homo sapiens",
        body_weight_kg=70.0,
        liver_weight_g=1800.0,
        hepatic_blood_flow_l_h=90.0,
        microsomal_protein_mg_per_g_liver=45.0,
        source="Test source",
    )
    clint_low = InVitroAdmeRecord(
        record_id="ADME-MIC-LOW",
        compound_id="PW-TEST",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=["Research use only. Not a dosing recommendation."],
        value=0.01,  # very low
        unit="uL/min/mg_protein",
    )
    pred_low = predict_hepatic_clearance_well_stirred(
        clint_low,
        physiology=phys,
        fu_plasma=0.01,
        fu_inc=1.0,
        blood_to_plasma_ratio=1.0,
    )
    assert pred_low.hepatic_extraction_ratio < 0.001
    # CLH,b is almost exactly fu_b * CLint,u,H
    assert pred_low.hepatic_blood_clearance_l_h == pytest.approx(
        pred_low.fu_blood * pred_low.clint_u_hepatic_l_h, rel=1e-3
    )

    # High extraction: CLH,b -> QH
    clint_high = InVitroAdmeRecord(
        record_id="ADME-MIC-HIGH",
        compound_id="PW-TEST",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=["Research use only. Not a dosing recommendation."],
        value=5000.0,  # very high
        unit="uL/min/mg_protein",
    )
    pred_high = predict_hepatic_clearance_well_stirred(
        clint_high,
        physiology=phys,
        fu_plasma=1.0,
        fu_inc=1.0,
        blood_to_plasma_ratio=1.0,
    )
    # Extraction ratio approaches 1.0
    assert 0.99 < pred_high.hepatic_extraction_ratio < 1.0
    assert pred_high.hepatic_blood_clearance_l_h == pytest.approx(phys.hepatic_blood_flow_l_h, rel=1e-2)


def test_fold_metrics_and_summary_diagnostics() -> None:
    """Verify fold metrics, AFE, AAFE, and outlier fractions."""
    # 1. Exact match
    m_exact = compute_fold_metrics(10.0, 10.0)
    assert m_exact.ratio == 1.0
    assert m_exact.fold_error == 1.0
    assert m_exact.log10_error == 0.0
    assert m_exact.within_2_fold is True
    assert m_exact.within_3_fold is True
    assert m_exact.outlier_2_fold is False
    assert m_exact.bias_direction == BiasDirection.NONE

    # 2. Over-prediction by 2-fold
    m_over = compute_fold_metrics(20.0, 10.0)
    assert m_over.ratio == 2.0
    assert m_over.fold_error == 2.0
    assert m_over.log10_error == pytest.approx(math.log10(2.0), rel=1e-9)
    assert m_over.within_2_fold is True
    assert m_over.bias_direction == BiasDirection.OVER

    # 3. Under-prediction by 2-fold (symmetric fold error is 2.0)
    m_under = compute_fold_metrics(5.0, 10.0)
    assert m_under.ratio == 0.5
    assert m_under.fold_error == 2.0
    assert m_under.log10_error == pytest.approx(-math.log10(2.0), rel=1e-9)
    assert m_under.within_2_fold is True
    assert m_under.bias_direction == BiasDirection.UNDER

    # 4. Outliers
    m_out = compute_fold_metrics(25.0, 10.0)  # 2.5-fold
    assert m_out.within_2_fold is False
    assert m_out.outlier_2_fold is True
    assert m_out.within_3_fold is True

    m_far = compute_fold_metrics(40.0, 10.0)  # 4-fold
    assert m_far.within_3_fold is False

    # Summary diagnostics: AFE vs AAFE
    # 1 entry 2-fold over, 1 entry 2-fold under
    # AFE = 10^( (log10(2) - log10(2)) / 2 ) = 10^0 = 1.0 (cancels bias)
    # AAFE = 10^( (|log10(2)| + |-log10(2)|) / 2 ) = 10^(log10(2)) = 2.0 (does not cancel)
    entries = [
        PredictionAuditEntry(
            prior_evidence_id="P1",
            estimate_evidence_id="E1",
            prior_source_project="theraweave",
            compound_id="C1",
            measurement_type="clearance",
            prior_claim_type="clearance",
            estimate_claim_type="CL",
            predicted_value_original=20.0,
            predicted_unit_original="L/h",
            predicted_value=20.0,
            estimated_value=10.0,
            unit="L/h",
            conversion="same_unit",
            ratio=m_over.ratio,
            fold_error=m_over.fold_error,
            log10_error=m_over.log10_error,
            within_2_fold=m_over.within_2_fold,
            within_3_fold=m_over.within_3_fold,
            outlier_2_fold=m_over.outlier_2_fold,
            bias_direction=m_over.bias_direction,
        ),
        PredictionAuditEntry(
            prior_evidence_id="P2",
            estimate_evidence_id="E2",
            prior_source_project="theraweave",
            compound_id="C1",
            measurement_type="clearance",
            prior_claim_type="clearance",
            estimate_claim_type="CL",
            predicted_value_original=5.0,
            predicted_unit_original="L/h",
            predicted_value=5.0,
            estimated_value=10.0,
            unit="L/h",
            conversion="same_unit",
            ratio=m_under.ratio,
            fold_error=m_under.fold_error,
            log10_error=m_under.log10_error,
            within_2_fold=m_under.within_2_fold,
            within_3_fold=m_under.within_3_fold,
            outlier_2_fold=m_under.outlier_2_fold,
            bias_direction=m_under.bias_direction,
        ),
    ]
    summaries = summarize_audit_entries(entries)
    assert len(summaries) == 1
    s = summaries[0]
    assert s.measurement_type == "clearance"
    assert s.n_comparisons == 2
    assert s.afe == pytest.approx(1.0, rel=1e-9)
    assert s.aafe == pytest.approx(2.0, rel=1e-9)
    assert s.fraction_within_2_fold == 1.0
    assert s.fraction_within_3_fold == 1.0
    assert s.n_outlier_2_fold == 0


def test_audit_study_predictions_end_to_end(repo_root) -> None:
    """Test full calibration audit with body-weight conversion, semantic isolation, and immutability."""
    study = build_study()
    nca_results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    # SUBJ-01: BW = 0.25 kg, CL in vivo estimate = 0.05 L/h (0.20 L/h/kg)
    # SUBJ-02: BW = 0.30 kg, CL in vivo estimate = 0.06 L/h (0.20 L/h/kg)
    estimates = []
    for row in nca_results:
        estimates.extend(evidence_from_nca(row, run_id="run-nca"))

    # Create priors:
    # 1. Prior for PW-0001 clearance in L/h/kg (0.22 L/h/kg) -> body-weight normalized
    #    For SUBJ-01 (0.25 kg): predicted = 0.22 * 0.25 = 0.055 L/h vs estimated 0.05 L/h -> ratio = 1.1
    #    For SUBJ-02 (0.30 kg): predicted = 0.22 * 0.30 = 0.066 L/h vs estimated 0.06 L/h -> ratio = 1.1
    prior_cl_bw = EvidenceRecord(
        evidence_id="PRIOR-CL-BW",
        source_project=ProjectName.THERAWEAVE,
        source_run_id="tw-run-1",
        entity_type="CandidateCompound",
        entity_id="PW-0001",
        measurement_type="clearance",
        evidence_type=EvidenceType.PREDICTED,
        qc_status=EvidenceQc.NOT_APPLICABLE,
        claim_type="clearance",
        numeric_value=0.22,
        unit="L/h/kg",
        limitations=["Predicted clearance"],
    )

    # 2. Prior for apparent_clearance (CL/F) -> must NOT match IV clearance estimates!
    prior_clf = EvidenceRecord(
        evidence_id="PRIOR-CLF",
        source_project=ProjectName.THERAWEAVE,
        source_run_id="tw-run-1",
        entity_type="CandidateCompound",
        entity_id="PW-0001",
        measurement_type="apparent_clearance",
        evidence_type=EvidenceType.PREDICTED,
        qc_status=EvidenceQc.NOT_APPLICABLE,
        claim_type="CL_over_F",
        numeric_value=0.08,
        unit="L/h",
        limitations=["Predicted oral clearance"],
    )

    # 3. Prior for unknown compound PW-9999 -> skipped
    prior_unknown = EvidenceRecord(
        evidence_id="PRIOR-UNKNOWN",
        source_project=ProjectName.THERAWEAVE,
        source_run_id="tw-run-1",
        entity_type="CandidateCompound",
        entity_id="PW-9999",
        measurement_type="clearance",
        evidence_type=EvidenceType.PREDICTED,
        qc_status=EvidenceQc.NOT_APPLICABLE,
        claim_type="clearance",
        numeric_value=1.0,
        unit="L/h",
        limitations=["Predicted clearance"],
    )

    priors = [prior_cl_bw, prior_clf, prior_unknown]

    report = audit_study_predictions(study, estimates, priors=priors, audit_id="TEST-AUDIT")
    assert report.audit_id == "TEST-AUDIT"
    assert report.n_priors == 3

    # CL/F and PW-9999 must be skipped
    skipped_ids = {s.prior_evidence_id for s in report.skipped}
    assert "PRIOR-CLF" in skipped_ids
    assert "PRIOR-UNKNOWN" in skipped_ids

    # PRIOR-CL-BW matched 2 subjects
    assert len(report.entries) == 2
    for e in report.entries:
        assert e.prior_evidence_id == "PRIOR-CL-BW"
        assert e.conversion == "body_weight"
        assert e.ratio == pytest.approx(1.1, rel=1e-4)
        assert e.within_2_fold is True

    # Summaries
    assert len(report.summaries) == 1
    assert report.summaries[0].measurement_type == "clearance"
    assert report.summaries[0].afe == pytest.approx(1.1, rel=1e-4)
    assert report.overall_aafe == pytest.approx(1.1, rel=1e-4)
    assert report.overall_fraction_within_2_fold == 1.0


def test_audit_contract_violations() -> None:
    """Non-predicted priors or non-estimated posteriors violate contracts."""
    valid_prior = EvidenceRecord(
        evidence_id="PRIOR-1",
        source_project=ProjectName.THERAWEAVE,
        source_run_id="tw-1",
        entity_type="CandidateCompound",
        entity_id="PW-0001",
        measurement_type="clearance",
        evidence_type=EvidenceType.PREDICTED,
        qc_status=EvidenceQc.NOT_APPLICABLE,
        claim_type="clearance",
        numeric_value=0.1,
        unit="L/h",
        limitations=["Research use only."],
    )
    valid_estimate = EvidenceRecord(
        evidence_id="EST-1",
        source_project=ProjectName.PHARMACOWEAVE,
        source_run_id="pw-1",
        entity_type="NcaParameterSet",
        entity_id="SUBJ-01:PW-0001",
        measurement_type="clearance",
        evidence_type=EvidenceType.ESTIMATED,
        qc_status=EvidenceQc.NOT_REVIEWED,
        claim_type="clearance",
        numeric_value=0.1,
        unit="L/h",
        derived_from=["nca:1"],
        limitations=["Research use only."],
    )

    # Passing observed evidence as prior
    observed_as_prior = valid_prior.model_copy(update={"evidence_type": EvidenceType.OBSERVED})
    with pytest.raises(ContractViolation, match="only evidence_type=predicted"):
        audit_evidence_predictions([observed_as_prior], [valid_estimate])

    # Passing predicted evidence as estimate
    predicted_as_estimate = valid_estimate.model_copy(update={"evidence_type": EvidenceType.PREDICTED})
    with pytest.raises(ContractViolation, match="must be evidence_type=estimated"):
        audit_evidence_predictions([valid_prior], [predicted_as_estimate])
