"""Cross-project evidence and identity records.

An evidence record cites a result from one project. It does not become a
concentration row. Identity links are explicit and carry a review status;
matching subject-id text is not a link.
"""

from enum import StrEnum

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Identifier, Record
from pharmacoweave.version import SCHEMA_VERSION


class EvidenceType(StrEnum):
    OBSERVED = "observed"
    ESTIMATED = "estimated"
    PREDICTED = "predicted"
    LITERATURE = "literature"
    HYPOTHESIS = "hypothesis"


class ReviewStatus(StrEnum):
    UNREVIEWED = "unreviewed"
    ACCEPTED = "accepted"
    REJECTED = "rejected"


class ProjectName(StrEnum):
    PHARMACOWEAVE = "pharmacoweave"
    OMICSWEAVE = "omicsweave"
    THERAWEAVE = "theraweave"
    EPHYSWEAVE = "ephysweave"
    EXTERNAL = "external"


class EvidenceQc(StrEnum):
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    NOT_REVIEWED = "not_reviewed"
    NOT_APPLICABLE = "not_applicable"


RESEARCH_USE_LIMITATION = "Research use only. Not a dosing recommendation."
TOX_RESEARCH_LIMITATION = "Research use only. Not a clinical dosing recommendation or a safety claim."


class EvidenceRecord(Record):
    """One claim, with its evidence grade kept on the record."""

    evidence_id: Identifier
    source_project: ProjectName
    source_run_id: str = Field(min_length=1)
    entity_type: str = Field(min_length=1)
    entity_id: str = Field(min_length=1)
    measurement_type: str = Field(min_length=1)
    evidence_type: EvidenceType
    qc_status: EvidenceQc
    claim_type: str = Field(min_length=1)
    schema_name: str = "pharmacology-evidence"
    schema_version: str = SCHEMA_VERSION
    source_artifact_id: str | None = None
    unit: str | None = None
    time_context: str | None = None
    applicability: str | None = None
    limitations: list[str] = Field(default_factory=list)
    derived_from: list[str] = Field(default_factory=list)
    subject_id: str | None = None
    session_id: str | None = None
    timebase_id: str | None = None
    sync_error_ms: float | None = None
    numeric_value: float | None = None

    @model_validator(mode="after")
    def contract_version(self) -> "EvidenceRecord":
        if self.schema_name != "pharmacology-evidence":
            raise ValueError("schema_name must be pharmacology-evidence.")
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        if self.sync_error_ms is not None and self.sync_error_ms < 0:
            raise ValueError("sync_error_ms cannot be negative.")
        return self


class EvidenceExchange(Record):
    """A file-level bundle of evidence records."""

    records: list[EvidenceRecord] = Field(default_factory=list)
    schema_name: str = "pharmacology-evidence-exchange"
    schema_version: str = SCHEMA_VERSION

    @model_validator(mode="after")
    def contract_version(self) -> "EvidenceExchange":
        if self.schema_name != "pharmacology-evidence-exchange":
            raise ValueError("schema_name must be pharmacology-evidence-exchange.")
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        ids = [record.evidence_id for record in self.records]
        if len(ids) != len(set(ids)):
            raise ValueError("evidence_id values must be unique inside an exchange.")
        return self


class IdentityLink(Record):
    """A reviewed correspondence between a foreign subject id and a study subject."""

    link_id: Identifier
    source_project: ProjectName
    source_study_id: str = Field(min_length=1)
    source_subject_id: str = Field(min_length=1)
    canonical_subject_id: Identifier
    match_method: str = Field(min_length=3)
    review_status: ReviewStatus
    notes: str | None = None


class SourceFile(Record):
    """Hash of one file that was read to build a study."""

    relative_path: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    role: str = Field(min_length=1)
