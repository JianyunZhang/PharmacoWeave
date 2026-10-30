"""In vitro ADME records stay in vitro."""

import pytest
from pydantic import ValidationError

from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.dmpk import (
    AdmeAssayKind,
    AdmeStatus,
    InVitroAdmeRecord,
    scale_to_in_vivo_clearance,
    validate_adme_records,
)
from pharmacoweave.errors import AnalysisError


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
