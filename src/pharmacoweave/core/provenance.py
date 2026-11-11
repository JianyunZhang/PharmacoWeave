"""Result-bundle manifest.

A bundle written without analyses records validated study data with
``analysis_executed=False`` and lists unexecuted analysis roles under
``omitted`` so a missing file is never read as an empty analysis.
"""

from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Identifier, Record, UtcDateTime
from pharmacoweave.version import SCHEMA_VERSION


class Producer(Record):
    """The software that wrote a bundle."""

    software: Literal["PharmacoWeave"] = "PharmacoWeave"
    version: str = Field(min_length=1)


class BundleFile(Record):
    """One file inside a bundle, hashed from the bytes that were written."""

    path: str = Field(min_length=1)
    role: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    media_type: str = Field(min_length=1)


class OmittedArtifact(Record):
    """A planned artifact that this run deliberately did not write."""

    role: str = Field(min_length=1)
    reason: str = Field(min_length=1)


class SourceBundleRef(Record):
    """A pointer to an upstream bundle. The upstream payload is not copied."""

    schema_name: str = Field(min_length=1)
    schema_version: str = Field(min_length=1)
    source_project: str = Field(min_length=1)
    source_run_id: str = Field(min_length=1)
    uri: str = Field(min_length=1)


class ResultBundleManifest(Record):
    """Index of a pharmacology result bundle."""

    run_id: Identifier
    study_id: Identifier
    created_at: UtcDateTime
    producer: Producer
    files: list[BundleFile]
    data: dict[str, str]
    omitted: list[OmittedArtifact]
    analysis_executed: bool
    source_bundles: list[SourceBundleRef] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    schema_name: Literal["pharmacology-result-bundle"] = "pharmacology-result-bundle"
    schema_version: Literal["0.1.0"] = "0.1.0"
    quality_status: Literal["research_only"] = "research_only"
    provenance: str = "provenance/source_files.json"

    @model_validator(mode="after")
    def analysis_flag_matches_files(self) -> "ResultBundleManifest":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        analysis_roles = {
            "nca_parameters",
            "model_parameters",
            "predictions",
            "residuals",
            "diagnostics",
            "emax_parameters",
            "emax_predictions",
            "derived_observations",
            "ivive_predictions",
            "prediction_audit",
            "exposure_simulation",
            "accepted_record",
        }
        listed = analysis_roles & set(self.data)
        if self.analysis_executed and not listed:
            raise ValueError("analysis_executed requires analysis artifacts to be listed in data.")
        if not self.analysis_executed and listed:
            raise ValueError("Analysis artifacts cannot be listed when analysis_executed is false.")
        return self
