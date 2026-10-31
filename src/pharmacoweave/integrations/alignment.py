"""Multi-scale electrophysiology and omics endpoint alignment with study exposures.

Enforces cross-project boundary contracts:
1. Subject mapping must have review_status='accepted' in the study's identity_map.
2. Clocks must be calibrated with explicit timebase_id matching study timebases.
3. Synchronization error (sync_error_ms) must not exceed the allowed threshold.
4. Omics metabolic flux is not a concentration and cannot be ingested as one.
5. Preserves end-to-end provenance: dose_event -> exposure -> foreign PD endpoint -> source artifact.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from pydantic import Field, model_validator

if TYPE_CHECKING:
    from pharmacoweave.core.study import PharmacologyStudy

from pharmacoweave.core.entities import Identifier, Record
from pharmacoweave.core.evidence import (
    RESEARCH_USE_LIMITATION,
    EvidenceExchange,
    EvidenceRecord,
    ProjectName,
    ReviewStatus,
)
from pharmacoweave.core.units import CONC_AMOUNT, CONC_MASS, UnitError, parse_unit
from pharmacoweave.errors import AlignmentError, ContractViolation
from pharmacoweave.version import SCHEMA_VERSION


class AlignedEndpointRecord(Record):
    """A cross-project PD/omics endpoint successfully aligned to a study exposure."""

    aligned_id: Identifier
    source_project: ProjectName
    source_entity_id: str = Field(min_length=1)
    source_artifact_id: str = Field(min_length=1)
    canonical_subject_id: Identifier
    dose_event_id: Identifier
    exposure_ref: str = Field(min_length=1)
    timebase_id: Identifier
    sync_error_ms: float = Field(ge=0.0)
    measurement_type: str = Field(min_length=1)
    numeric_value: float
    unit: str = Field(min_length=1)
    traceability_path: str = Field(min_length=1)
    derived_from: list[str] = Field(min_length=1)
    limitations: list[str] = Field(min_length=1)
    evidence_grade: Literal["observed"] = "observed"
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def validate_aligned_record(self) -> AlignedEndpointRecord:
        if RESEARCH_USE_LIMITATION not in self.limitations:
            raise ValueError(f"limitations must include {RESEARCH_USE_LIMITATION!r}.")
        return self


class RejectedAlignment(Record):
    """An alignment attempt that was rejected by boundary or synchrony rules."""

    source_project: ProjectName
    source_entity_id: str = Field(min_length=1)
    reason: str = Field(min_length=1)
    sync_error_ms: float | None = None
    subject_id: str | None = None
    timebase_id: str | None = None


class MultiScaleAlignmentReport(Record):
    """Full report of cross-project multi-scale endpoint alignment."""

    study_id: Identifier
    max_allowed_sync_error_ms: float = Field(gt=0.0)
    aligned_endpoints: list[AlignedEndpointRecord] = Field(default_factory=list)
    rejected_alignments: list[RejectedAlignment] = Field(default_factory=list)
    schema_name: Literal["pharmacology-multiscale-alignment"] = "pharmacology-multiscale-alignment"
    schema_version: str = SCHEMA_VERSION
    research_use_only: Literal[True] = True


def _check_flux_boundary(record: EvidenceRecord) -> None:
    """Refuse attempts to classify metabolic flux as a pharmacokinetic concentration."""
    measurement = record.measurement_type.lower()
    claim = record.claim_type.lower()
    entity = record.entity_type.lower()

    is_flux = "flux" in measurement or "flux" in claim or "flux" in entity
    if is_flux:
        if record.unit is not None:
            try:
                dims = parse_unit(record.unit).dimensions
                if dims in {CONC_MASS, CONC_AMOUNT}:
                    raise ContractViolation(
                        f"Metabolic flux record {record.evidence_id} cannot be represented "
                        f"with concentration dimensions ({record.unit}). Flux is not concentration."
                    )
            except UnitError:
                pass
        if claim in {"concentration", "pk_concentration"} or entity in {"concentration", "analyte"}:
            raise ContractViolation(
                f"Metabolic flux record {record.evidence_id} cannot be claimed as "
                f"a concentration observation. Flux is not concentration."
            )


def align_cross_project_endpoint(
    *,
    study: PharmacologyStudy,
    foreign_record: EvidenceRecord,
    dose_event_id: str,
    exposure_ref: str,
    max_allowed_sync_error_ms: float = 1000.0,
) -> AlignedEndpointRecord:
    """Align a single EphysWeave or OmicsWeave endpoint to a study dose and exposure.

    Raises:
        ContractViolation: if an Omics flux is claimed as a concentration.
        AlignmentError: if identity link is missing/unaccepted, clocks do not match,
                        or sync_error_ms exceeds the allowed threshold.
    """
    _check_flux_boundary(foreign_record)

    if not foreign_record.subject_id:
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} lacks subject_id; cannot align."
        )

    # 1. Identity link check
    matching_link = next(
        (
            link
            for link in study.identity_links
            if link.source_project == foreign_record.source_project
            and link.source_subject_id == foreign_record.subject_id
        ),
        None,
    )
    if matching_link is None:
        raise AlignmentError(
            f"No identity link found for {foreign_record.source_project.value} "
            f"subject {foreign_record.subject_id!r} in study identity_map."
        )
    if matching_link.review_status != ReviewStatus.ACCEPTED:
        raise AlignmentError(
            f"Identity link {matching_link.link_id} for subject {foreign_record.subject_id!r} "
            f"has review_status={matching_link.review_status.value!r}; "
            "cross-project alignment requires review_status='accepted'."
        )
    canonical_subject_id = matching_link.canonical_subject_id

    # 2. Dose event and subject consistency
    dose_event = next(
        (d for d in study.dose_events if d.dose_event_id == dose_event_id),
        None,
    )
    if dose_event is None:
        raise AlignmentError(
            f"Referenced dose_event_id {dose_event_id!r} does not exist in study dose_events."
        )
    if dose_event.subject_id != canonical_subject_id:
        raise AlignmentError(
            f"Dose event {dose_event_id} was administered to subject {dose_event.subject_id}, "
            f"which does not match canonical subject {canonical_subject_id}."
        )

    # 3. Timebase verification
    study_timebase_ids = {tb.timebase_id for tb in study.metadata.timebases}
    if not foreign_record.timebase_id:
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} does not declare timebase_id."
        )
    if foreign_record.timebase_id not in study_timebase_ids:
        raise AlignmentError(
            f"Timebase {foreign_record.timebase_id!r} declared by foreign record "
            f"is not registered in study timebases: {sorted(study_timebase_ids)}."
        )

    # 4. Clock synchronization error verification
    if foreign_record.sync_error_ms is None:
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} lacks sync_error_ms. "
            "Cross-scale alignment requires explicit clock synchrony estimation."
        )
    if foreign_record.sync_error_ms < 0:
        raise AlignmentError(
            f"sync_error_ms cannot be negative ({foreign_record.sync_error_ms} ms)."
        )
    if foreign_record.sync_error_ms > max_allowed_sync_error_ms:
        raise AlignmentError(
            f"Clock synchronization error {foreign_record.sync_error_ms:.2f} ms exceeds "
            f"maximum allowed threshold ({max_allowed_sync_error_ms:.2f} ms). Alignment refused."
        )

    # 5. Source artifact presence
    if not foreign_record.source_artifact_id or not foreign_record.source_artifact_id.strip():
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} lacks source_artifact_id. "
            "End-to-end traceability requires an upstream artifact hash or URI."
        )

    if foreign_record.numeric_value is None:
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} lacks numeric_value; cannot align."
        )
    if not foreign_record.unit:
        raise AlignmentError(
            f"Evidence record {foreign_record.evidence_id} lacks unit; cannot align."
        )

    # 6. Build full provenance chain: dose_event -> exposure -> foreign PD endpoint -> source artifact
    trace_path = (
        f"dose:{dose_event_id} -> exposure:{exposure_ref} -> "
        f"{foreign_record.source_project.value}:{foreign_record.entity_id} -> "
        f"artifact:{foreign_record.source_artifact_id}"
    )
    derived_from = [
        f"dose:{dose_event_id}",
        f"exposure:{exposure_ref}",
        f"{foreign_record.source_project.value}:{foreign_record.entity_id}",
        f"artifact:{foreign_record.source_artifact_id}",
    ]
    if foreign_record.derived_from:
        derived_from.extend(foreign_record.derived_from)

    aligned_id = f"aligned_{foreign_record.source_project.value}_{foreign_record.evidence_id}"

    return AlignedEndpointRecord(
        aligned_id=aligned_id,
        source_project=foreign_record.source_project,
        source_entity_id=foreign_record.entity_id,
        source_artifact_id=foreign_record.source_artifact_id,
        canonical_subject_id=canonical_subject_id,
        dose_event_id=dose_event_id,
        exposure_ref=exposure_ref,
        timebase_id=foreign_record.timebase_id,
        sync_error_ms=foreign_record.sync_error_ms,
        measurement_type=foreign_record.measurement_type,
        numeric_value=foreign_record.numeric_value,
        unit=foreign_record.unit,
        traceability_path=trace_path,
        derived_from=derived_from,
        limitations=[RESEARCH_USE_LIMITATION, *foreign_record.limitations],
    )


def align_study_cross_project_evidence(
    *,
    study: PharmacologyStudy,
    exchange: EvidenceExchange,
    dose_exposure_mapping: dict[str, tuple[str, str]],
    max_allowed_sync_error_ms: float = 1000.0,
) -> MultiScaleAlignmentReport:
    """Align a collection of foreign evidence records with study exposures.

    dose_exposure_mapping maps foreign record evidence_id -> (dose_event_id, exposure_ref).
    Records failing alignment are captured under rejected_alignments with full diagnostic reasons.
    """
    aligned: list[AlignedEndpointRecord] = []
    rejected: list[RejectedAlignment] = []

    for record in exchange.records:
        if record.source_project not in {ProjectName.EPHYSWEAVE, ProjectName.OMICSWEAVE}:
            continue

        mapping = dose_exposure_mapping.get(record.evidence_id)
        if mapping is None:
            rejected.append(
                RejectedAlignment(
                    source_project=record.source_project,
                    source_entity_id=record.entity_id,
                    reason=f"No dose/exposure mapping provided for evidence {record.evidence_id}.",
                    sync_error_ms=record.sync_error_ms,
                    subject_id=record.subject_id,
                    timebase_id=record.timebase_id,
                )
            )
            continue

        dose_event_id, exposure_ref = mapping
        try:
            item = align_cross_project_endpoint(
                study=study,
                foreign_record=record,
                dose_event_id=dose_event_id,
                exposure_ref=exposure_ref,
                max_allowed_sync_error_ms=max_allowed_sync_error_ms,
            )
            aligned.append(item)
        except (AlignmentError, ContractViolation) as exc:
            rejected.append(
                RejectedAlignment(
                    source_project=record.source_project,
                    source_entity_id=record.entity_id,
                    reason=str(exc),
                    sync_error_ms=record.sync_error_ms,
                    subject_id=record.subject_id,
                    timebase_id=record.timebase_id,
                )
            )

    return MultiScaleAlignmentReport(
        study_id=study.metadata.study_id,
        max_allowed_sync_error_ms=max_allowed_sync_error_ms,
        aligned_endpoints=aligned,
        rejected_alignments=rejected,
    )
