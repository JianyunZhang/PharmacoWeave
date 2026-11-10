"""Rules that keep evidence grades from being promoted across projects."""

import json
from dataclasses import dataclass
from pathlib import Path

from pharmacoweave.core.evidence import EvidenceExchange, EvidenceQc, EvidenceRecord, EvidenceType, ProjectName
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import ContractViolation, IngestError

_PK_MEASUREMENTS = frozenset(
    {
        "concentration",
        "clearance",
        "apparent_clearance",
        "auc",
        "half_life",
        "volume",
        "apparent_volume",
        "ka",
        "renal_clearance",
        "amount_excreted",
    }
)
_EPHYS_ENDPOINTS = frozenset({"neural_pd_endpoint", "behavioral_endpoint"})
_NON_MEASUREMENT = frozenset({EvidenceType.PREDICTED, EvidenceType.LITERATURE, EvidenceType.HYPOTHESIS})


@dataclass(frozen=True)
class EphysEndpointRef:
    """A timed EphysWeave endpoint that has not been aligned to a dose clock."""

    evidence_id: str
    measurement_type: str
    subject_id: str
    session_id: str
    timebase_id: str
    sync_error_ms: float
    value: float
    unit: str


def screen_record(record: EvidenceRecord) -> None:
    """Raise when a record would smuggle PK, predictions, or unsynced endpoints."""

    if record.evidence_type is EvidenceType.OBSERVED and record.measurement_type in _PK_MEASUREMENTS:
        raise ContractViolation(
            f"{record.evidence_id}: observed {record.measurement_type} belongs in the study "
            "observation table, not in an evidence exchange."
        )
    if record.source_project is ProjectName.OMICSWEAVE and record.measurement_type in _PK_MEASUREMENTS:
        raise ContractViolation(
            f"{record.evidence_id}: OmicsWeave {record.measurement_type} cannot be used as a PK measurement."
        )
    if (
        record.source_project is ProjectName.THERAWEAVE
        and record.evidence_type in {EvidenceType.OBSERVED, EvidenceType.ESTIMATED}
        and record.measurement_type in _PK_MEASUREMENTS
    ):
        raise ContractViolation(
            f"{record.evidence_id}: TheraWeave cannot supply observed or estimated PK. "
            "A predicted ADMET value has to stay evidence_type=predicted."
        )
    if record.evidence_type in _NON_MEASUREMENT and record.qc_status is EvidenceQc.ACCEPTED:
        raise ContractViolation(
            f"{record.evidence_id}: qc_status accepted is reserved for reviewed measurements. "
            "Predictions, literature, and hypotheses use not_applicable."
        )
    if record.source_project is ProjectName.EPHYSWEAVE and record.measurement_type in _EPHYS_ENDPOINTS:
        read_ephys_endpoint(record)
    if record.source_project is ProjectName.PHARMACOWEAVE and (
        not record.derived_from or any(not item.strip() for item in record.derived_from)
    ):
        raise ContractViolation(
            f"{record.evidence_id}: PharmacoWeave evidence requires non-empty derived_from links "
            "to the exposure or analysis records it came from."
        )


def screen_exchange(exchange: EvidenceExchange) -> None:
    """Screen every record in an exchange."""

    for record in exchange.records:
        screen_record(record)


def read_ephys_endpoint(record: EvidenceRecord) -> EphysEndpointRef:
    """Read an EphysWeave endpoint without shifting its clock or interpolating it."""

    if record.source_project is not ProjectName.EPHYSWEAVE:
        raise ContractViolation(f"{record.evidence_id} is not an EphysWeave record.")
    if record.measurement_type not in _EPHYS_ENDPOINTS:
        raise ContractViolation(
            f"{record.evidence_id} measurement_type must be neural_pd_endpoint or behavioral_endpoint."
        )
    if record.evidence_type is not EvidenceType.OBSERVED:
        raise ContractViolation(f"{record.evidence_id}: an EphysWeave endpoint must be an observed measurement.")
    required = {
        "subject_id": record.subject_id,
        "session_id": record.session_id,
        "timebase_id": record.timebase_id,
        "sync_error_ms": record.sync_error_ms,
        "numeric_value": record.numeric_value,
        "unit": record.unit,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        joined = ", ".join(missing)
        raise ContractViolation(
            f"{record.evidence_id} is missing {joined}. Clock alignment is not attempted when "
            "sync_error_ms or timebase_id is absent."
        )
    assert record.subject_id is not None
    assert record.session_id is not None
    assert record.timebase_id is not None
    assert record.sync_error_ms is not None
    assert record.numeric_value is not None
    assert record.unit is not None
    return EphysEndpointRef(
        evidence_id=record.evidence_id,
        measurement_type=record.measurement_type,
        subject_id=record.subject_id,
        session_id=record.session_id,
        timebase_id=record.timebase_id,
        sync_error_ms=record.sync_error_ms,
        value=record.numeric_value,
        unit=record.unit,
    )


def promote_to_observation(record: EvidenceRecord) -> None:
    """Refuse to turn evidence into an observation.

    Measured concentrations are loaded with the study tables. EphysWeave
    endpoints stay endpoint references until a later version aligns them to
    dose events under an explicit sync check.
    """

    raise ContractViolation(
        f"{record.evidence_id} cannot be promoted to an Observation. "
        "Evidence grade and project of origin stay on the evidence record."
    )


def attach_evidence_file(study: PharmacologyStudy, path: str | Path) -> PharmacologyStudy:
    """Attach a JSON evidence exchange after screening it.

    The observation table and compound table are left unchanged.
    """

    exchange = load_exchange(path)
    screen_exchange(exchange)
    return study.model_copy(update={"external_evidence": [*study.external_evidence, *exchange.records]})


def load_exchange(path: str | Path) -> EvidenceExchange:
    """Load one JSON evidence exchange from a file or a single-file directory."""

    selected = Path(path)
    if selected.is_dir():
        files = sorted(item for item in selected.glob("*.json") if item.is_file())
        if len(files) != 1:
            raise IngestError(f"{selected} must contain exactly one JSON evidence exchange, found {len(files)}.")
        selected = files[0]
    try:
        payload = json.loads(selected.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise IngestError(f"{selected} is not valid JSON: {exc}") from exc
    try:
        return EvidenceExchange.model_validate(payload)
    except Exception as exc:
        raise IngestError(f"{selected} does not match the evidence exchange contract: {exc}") from exc
