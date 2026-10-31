"""In vitro ADME records, binding assays, and exploratory IVIVE."""

import pytest
from pydantic import ValidationError

from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.dmpk import (
    AdmeAssayKind,
    AdmeStatus,
    InVitroAdmeRecord,
    InVitroAdmeSet,
    PhysiologicalScalingParams,
    human_starting_dose_from_ivive,
    predict_hepatic_clearance_well_stirred,
    scale_to_in_vivo_clearance,
    validate_adme_records,
)
from pharmacoweave.errors import AnalysisError, ContractViolation


def _record(**overrides: object) -> InVitroAdmeRecord:
    payload: dict[str, object] = {
        "record_id": "ADME-1",
        "compound_id": "PW-0001",
        "assay_kind": AdmeAssayKind.PERMEABILITY,
        "system": "Caco-2",
        "status": AdmeStatus.MEASURED,
        "limitations": [RESEARCH_USE_LIMITATION],
        "value": 1.0e-6,
        "unit": "cm/s",
    }
    payload.update(overrides)
    return InVitroAdmeRecord(**payload)  # type: ignore[arg-type]


def test_permeability_microsomal_and_hepatocyte_units() -> None:
    assert _record().unit == "cm/s"
    microsomal = _record(
        record_id="ADME-MIC",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        value=12.0,
        unit="uL/min/mg_protein",
    )
    hepatocyte = _record(
        record_id="ADME-HEP",
        assay_kind=AdmeAssayKind.HEPATOCYTE_CLEARANCE,
        system="human hepatocytes",
        value=4.0,
        unit="mL/min/million_cells",
    )
    validate_adme_records([_record(), microsomal, hepatocyte])
    assert microsomal.evidence_grade == "observed"


def test_plasma_clearance_and_missing_results_are_rejected() -> None:
    with pytest.raises(ValidationError, match="plasma clearance"):
        _record(
            assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
            unit="L/h",
            value=1.0,
        )
    with pytest.raises(ValidationError, match="uL/min/mg_protein"):
        _record(
            assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
            unit="uL/min/mg",
            value=1.0,
        )
    with pytest.raises(ValidationError, match="not carry a numeric"):
        _record(status=AdmeStatus.NOT_MEASURED, value=1.0, unit="cm/s")
    with pytest.raises(ValidationError, match="qc_reason"):
        _record(status=AdmeStatus.FAILED_QC, value=None, unit=None)
    failed = _record(status=AdmeStatus.FAILED_QC, value=None, unit=None, qc_reason="replicate disagreed")
    assert failed.value is None
    with pytest.raises(ValueError, match="unique"):
        validate_adme_records([_record(), _record()])


def test_in_vivo_scaling_is_not_performed() -> None:
    record = _record()
    with pytest.raises(AnalysisError, match="not scaled"):
        scale_to_in_vivo_clearance(record)


def test_fu_plasma_and_fu_mic_bounds_and_units() -> None:
    fu_p = _record(
        record_id="ADME-FUP",
        assay_kind=AdmeAssayKind.FU_PLASMA,
        system="human plasma",
        value=0.05,
        unit="fraction",
    )
    assert fu_p.canonical_value() == (0.05, "fraction")

    fu_pct = _record(
        record_id="ADME-FUP-PCT",
        assay_kind=AdmeAssayKind.FU_PLASMA,
        system="human plasma",
        value=5.0,
        unit="percent",
    )
    assert fu_pct.canonical_value() == (0.05, "fraction")

    # Greater than 1.0 fraction is rejected
    with pytest.raises(ValidationError, match="unbound fraction in"):
        _record(
            record_id="ADME-FUP-HIGH",
            assay_kind=AdmeAssayKind.FU_PLASMA,
            system="human plasma",
            value=1.5,
            unit="fraction",
        )

    # Greater than 100% is rejected
    with pytest.raises(ValidationError, match="unbound fraction in"):
        _record(
            record_id="ADME-FUP-HIGH-PCT",
            assay_kind=AdmeAssayKind.FU_PLASMA,
            system="human plasma",
            value=105.0,
            unit="percent",
        )

    # Zero or negative is rejected
    with pytest.raises(ValidationError, match="must be positive"):
        _record(
            record_id="ADME-FUP-ZERO",
            assay_kind=AdmeAssayKind.FU_PLASMA,
            system="human plasma",
            value=0.0,
            unit="fraction",
        )


def test_blood_to_plasma_ratio_units_and_bounds() -> None:
    bp_ratio = _record(
        record_id="ADME-BP",
        assay_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO,
        system="human whole blood",
        value=0.85,
        unit="ratio",
    )
    assert bp_ratio.canonical_value() == (0.85, "ratio")

    bp_fold = _record(
        record_id="ADME-BP-FOLD",
        assay_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO,
        system="human whole blood",
        value=1.2,
        unit="fold",
    )
    assert bp_fold.canonical_value() == (1.2, "ratio")

    with pytest.raises(ValidationError, match="must be one of"):
        _record(
            record_id="ADME-BP-BAD-UNIT",
            assay_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO,
            system="human whole blood",
            value=1.2,
            unit="mg/L",
        )


def test_adme_set_validation() -> None:
    rec1 = _record(record_id="ADME-1")
    rec2 = _record(record_id="ADME-2")
    adme_set = InVitroAdmeSet(records=[rec1, rec2])
    assert len(adme_set.records) == 2

    # Duplicate records rejected
    with pytest.raises(ValidationError, match="unique"):
        InVitroAdmeSet(records=[rec1, rec1])

    # Physiology profile
    phys = PhysiologicalScalingParams(
        profile_id="PHYS-HUMAN",
        species="Homo sapiens",
        body_weight_kg=70.0,
        liver_weight_g=1800.0,
        hepatic_blood_flow_l_h=90.0,
        microsomal_protein_mg_per_g_liver=45.0,
        source="Davies & Morris (1993)",
    )
    assert phys.profile_id == "PHYS-HUMAN"

    # Duplicate physiology profiles rejected
    with pytest.raises(ValidationError, match="profile_id"):
        InVitroAdmeSet(records=[rec1], physiology_profiles=[phys, phys])


def test_ivive_microsomal_and_hepatocyte_scaling() -> None:
    clint_mic = _record(
        record_id="ADME-MIC-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        value=20.0,
        unit="uL/min/mg_protein",
    )
    phys_human = PhysiologicalScalingParams(
        profile_id="PHYS-HUMAN-1",
        species="Homo sapiens",
        body_weight_kg=70.0,
        liver_weight_g=1800.0,
        hepatic_blood_flow_l_h=90.0,
        microsomal_protein_mg_per_g_liver=45.0,
        hepatocellularity_million_cells_per_g_liver=120.0,
        source="Davies & Morris (1993)",
    )
    fu_p = _record(
        record_id="ADME-FUP-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.FU_PLASMA,
        system="human plasma",
        value=0.1,
        unit="fraction",
    )
    fu_mic = _record(
        record_id="ADME-FUMIC-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.FU_MIC,
        system="human liver microsomes",
        value=0.8,
        unit="fraction",
    )
    bp = _record(
        record_id="ADME-BP-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.BLOOD_TO_PLASMA_RATIO,
        system="human whole blood",
        value=1.0,
        unit="ratio",
    )

    pred = predict_hepatic_clearance_well_stirred(
        clint_mic,
        physiology=phys_human,
        fu_plasma=fu_p,
        fu_inc=fu_mic,
        blood_to_plasma_ratio=bp,
    )
    assert pred.evidence_grade == "predicted"
    assert pred.compound_id == "PW-0001"
    assert pred.species == "Homo sapiens"
    assert pred.hepatic_plasma_clearance_l_h > 0
    assert 0 < pred.hepatic_extraction_ratio < 1
    assert any("Well-stirred" in a for a in pred.assumptions)
    assert any("Exploratory IVIVE" in lim for lim in pred.limitations)

    # Hepatocyte scaling with explicit scalar fu_inc
    clint_hep = _record(
        record_id="ADME-HEP-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.HEPATOCYTE_CLEARANCE,
        system="human hepatocytes",
        value=5.0,
        unit="uL/min/million_cells",
    )
    pred_hep = predict_hepatic_clearance_well_stirred(
        clint_hep,
        physiology=phys_human,
        fu_plasma=fu_p,
        fu_inc=1.0,
        blood_to_plasma_ratio=bp,
    )
    assert pred_hep.assay_kind == "hepatocyte_clearance"
    assert pred_hep.scaling_factor_unit == "million_cells/g_liver"

    # Human starting dose refusal
    with pytest.raises(ContractViolation, match="starting dose"):
        human_starting_dose_from_ivive(pred)


def test_ivive_guardrails_and_rejections() -> None:
    phys_no_mppgl = PhysiologicalScalingParams(
        profile_id="PHYS-NO-MIC",
        species="Homo sapiens",
        body_weight_kg=70.0,
        liver_weight_g=1800.0,
        hepatic_blood_flow_l_h=90.0,
        hepatocellularity_million_cells_per_g_liver=120.0,
        source="Test source",
    )
    clint_mic = _record(
        record_id="ADME-MIC-1",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="human liver microsomes",
        value=20.0,
        unit="uL/min/mg_protein",
    )
    with pytest.raises(AnalysisError, match="microsomal_protein_mg_per_g_liver"):
        predict_hepatic_clearance_well_stirred(
            clint_mic,
            physiology=phys_no_mppgl,
            fu_plasma=0.1,
            fu_inc=1.0,
            blood_to_plasma_ratio=1.0,
        )

    # Calling with permeability assay
    perm = _record(record_id="ADME-PERM", assay_kind=AdmeAssayKind.PERMEABILITY)
    with pytest.raises(AnalysisError, match="IVIVE needs microsomal_clearance or hepatocyte_clearance"):
        predict_hepatic_clearance_well_stirred(
            perm,
            physiology=phys_no_mppgl,
            fu_plasma=0.1,
            fu_inc=1.0,
            blood_to_plasma_ratio=1.0,
        )

    # Cross-compound binding borrowing rejected
    phys_with_mppgl = PhysiologicalScalingParams(
        profile_id="PHYS-WITH-MIC",
        species="Homo sapiens",
        body_weight_kg=70.0,
        liver_weight_g=1800.0,
        hepatic_blood_flow_l_h=90.0,
        microsomal_protein_mg_per_g_liver=45.0,
        source="Test source",
    )
    fu_p_diff = _record(
        record_id="ADME-FUP-DIFF",
        compound_id="PW-0002",
        assay_kind=AdmeAssayKind.FU_PLASMA,
        system="human plasma",
        value=0.1,
        unit="fraction",
    )
    with pytest.raises(AnalysisError, match="not borrowed across compounds"):
        predict_hepatic_clearance_well_stirred(
            clint_mic,
            physiology=phys_with_mppgl,
            fu_plasma=fu_p_diff,
            fu_inc=1.0,
            blood_to_plasma_ratio=1.0,
        )
