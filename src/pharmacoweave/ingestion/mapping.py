"""Declarative ingest mapping for heterogeneous laboratory tabular files.

Allows labs to map arbitrary column headers, delimiter styles, and format
peculiarities to canonical PharmacoWeave study tables while preserving the
cryptographic SHA-256 digests of all raw input files in study provenance.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any, Literal

import pyarrow.parquet as pq
import yaml
from pydantic import Field

from pharmacoweave.core.entities import (
    Analyte,
    Assay,
    Biosample,
    Compound,
    DoseEvent,
    Observation,
    Record,
    StudyMetadata,
    Subject,
)
from pharmacoweave.core.evidence import IdentityLink, SourceFile
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.errors import IngestError, ValidationFailed
from pharmacoweave.version import SCHEMA_VERSION

_MODELS: dict[str, type[Record]] = {
    "subjects": Subject,
    "compounds": Compound,
    "analytes": Analyte,
    "assays": Assay,
    "dosing_events": DoseEvent,
    "biosamples": Biosample,
    "observations": Observation,
    "identity_map": IdentityLink,
}

_NUMERIC_FLOAT_FIELDS = {
    "body_weight",
    "dose_value",
    "infusion_duration",
    "body_weight_normalized_value",
    "time_since_dose",
    "numeric_value",
    "loq_lower",
    "loq_upper",
    "molecular_weight",
}

_NUMERIC_INT_FIELDS = {
    "ncbi_taxon_id",
}

_BOOL_FIELDS = {
    "is_predose",
}


class TableMappingConfig(Record):
    """Declarative column and constant mapping for one study table."""

    source_file: str = Field(min_length=1)
    column_mapping: dict[str, str] = Field(default_factory=dict)
    constants: dict[str, Any] = Field(default_factory=dict)
    value_transforms: dict[str, dict[str, Any]] = Field(default_factory=dict)
    delimiter: str | None = None


class IngestMappingConfig(Record):
    """Configuration mapping raw lab files to a canonical PharmacologyStudy."""

    schema_name: Literal["pharmacology-ingest-mapping"] = "pharmacology-ingest-mapping"
    schema_version: str = SCHEMA_VERSION
    study_metadata: StudyMetadata
    tables: dict[str, TableMappingConfig]
    research_use_only: Literal[True] = True


def _compute_sha256(path: Path) -> str:
    """Calculate the hex SHA-256 digest of a file."""
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def _read_raw_rows(path: Path, delimiter: str | None = None) -> list[dict[str, Any]]:
    """Read raw tabular rows from a CSV, TSV, or Parquet file."""
    if not path.is_file():
        raise IngestError(f"Raw source file {path} does not exist.")

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        table = pq.read_table(path)
        return table.to_pylist()

    # Plain text CSV / TSV
    delim = delimiter
    if delim is None:
        delim = "\t" if suffix == ".tsv" else ","

    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        return list(reader)


def _cast_field(field_name: str, value: Any) -> Any:
    """Cast string values to appropriate target types if required."""
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "" or stripped.lower() in {"na", "null", "none", "nan"}:
            return None
        if field_name in _NUMERIC_FLOAT_FIELDS:
            try:
                return float(stripped)
            except ValueError as exc:
                raise IngestError(f"Could not parse float for field {field_name!r}: {stripped!r}") from exc
        if field_name in _NUMERIC_INT_FIELDS:
            try:
                return int(stripped)
            except ValueError as exc:
                raise IngestError(f"Could not parse int for field {field_name!r}: {stripped!r}") from exc
        if field_name in _BOOL_FIELDS:
            lower = stripped.lower()
            if lower in {"true", "1", "yes", "t", "y"}:
                return True
            if lower in {"false", "0", "no", "f", "n"}:
                return False
            raise IngestError(f"Could not parse boolean for field {field_name!r}: {stripped!r}")
        if field_name == "qc_flags":
            return [f.strip() for f in stripped.split(";") if f.strip()]
        return stripped
    return value


def load_ingest_mapping(mapping_path: str | Path) -> IngestMappingConfig:
    """Load and validate an ingest_mapping.yaml configuration."""
    path = Path(mapping_path)
    if not path.is_file():
        raise IngestError(f"Mapping configuration file {path} not found.")

    raw_yaml = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw_yaml, dict):
        raise IngestError(f"Invalid YAML structure in {path}: expected a dictionary.")

    try:
        return IngestMappingConfig.model_validate(raw_yaml)
    except Exception as exc:
        raise IngestError(f"Failed to validate ingest mapping schema in {path}: {exc}") from exc


def ingest_study_from_mapping(
    mapping: IngestMappingConfig | str | Path,
    *,
    base_dir: str | Path | None = None,
    mode: str = "research",
    validate_schema: bool = True,
) -> PharmacologyStudy:
    """Construct and validate a PharmacologyStudy from raw files using declarative mapping.

    The SHA-256 digest of every referenced raw laboratory file is recorded in
    study.source_files for reproducible provenance tracking.
    """
    if isinstance(mapping, (str, Path)):
        mapping_path = Path(mapping)
        config = load_ingest_mapping(mapping_path)
        actual_base_dir = Path(base_dir) if base_dir is not None else mapping_path.parent
    else:
        config = mapping
        if base_dir is None:
            raise IngestError("base_dir is required when passing an IngestMappingConfig object directly.")
        actual_base_dir = Path(base_dir)

    source_files: list[SourceFile] = []
    seen_files: dict[str, str] = {}  # rel_path -> sha256
    loaded_entities: dict[str, list[Record]] = {table_name: [] for table_name in _MODELS}

    for table_name, table_config in config.tables.items():
        if table_name not in _MODELS:
            raise IngestError(
                f"Unknown study table {table_name!r}. Allowed tables: {sorted(_MODELS.keys())}"
            )
        model_cls = _MODELS[table_name]
        raw_file_path = actual_base_dir / table_config.source_file
        if not raw_file_path.is_file():
            raise IngestError(
                f"Source file {table_config.source_file!r} for table {table_name!r} "
                f"not found at {raw_file_path}."
            )

        # Track cryptographic hash
        rel_str = table_config.source_file
        if rel_str not in seen_files:
            file_hash = _compute_sha256(raw_file_path)
            seen_files[rel_str] = file_hash
            source_files.append(
                SourceFile(
                    relative_path=rel_str,
                    sha256=file_hash,
                    role=f"raw_{table_name}",
                )
            )

        raw_rows = _read_raw_rows(raw_file_path, delimiter=table_config.delimiter)
        for row_idx, raw_row in enumerate(raw_rows):
            entity_data: dict[str, Any] = {}

            # 1. Apply column remapping
            for raw_col, target_field in table_config.column_mapping.items():
                if raw_col in raw_row:
                    val = raw_row[raw_col]
                    # Apply value transforms if specified
                    if target_field in table_config.value_transforms:
                        transform_map = table_config.value_transforms[target_field]
                        val = transform_map.get(str(val), val)
                    entity_data[target_field] = _cast_field(target_field, val)

            # 2. Apply constants for any missing target fields
            for target_field, const_val in table_config.constants.items():
                if target_field not in entity_data or entity_data[target_field] is None:
                    entity_data[target_field] = _cast_field(target_field, const_val)

            # Entity-specific adaptations
            if table_name == "subjects":
                if "body_weight" in entity_data and isinstance(entity_data["body_weight"], (int, float)):
                    bw_unit = entity_data.pop("body_weight_unit", "kg")
                    entity_data["body_weight"] = {"value": entity_data.pop("body_weight"), "unit": bw_unit}
                elif "body_weight_value" in entity_data:
                    bw_val = entity_data.pop("body_weight_value")
                    bw_unit = entity_data.pop("body_weight_unit", "kg")
                    entity_data["body_weight"] = None if bw_val is None else {"value": bw_val, "unit": bw_unit}

            if table_name == "observations" and "observation_id" not in entity_data:
                sample_part = entity_data.get("sample_id", "obs")
                entity_data["observation_id"] = f"{sample_part}_row{row_idx + 1}"

            # Deduplicate dimension tables if multiple rows declare the same entity
            if table_name == "compounds":
                cid = entity_data.get("compound_id")
                if any(getattr(c, "compound_id", None) == cid for c in loaded_entities["compounds"]):
                    continue
            if table_name == "analytes":
                aid = entity_data.get("analyte_id")
                if any(getattr(a, "analyte_id", None) == aid for a in loaded_entities["analytes"]):
                    continue
            if table_name == "assays":
                asid = entity_data.get("assay_id")
                if any(getattr(a, "assay_id", None) == asid for a in loaded_entities["assays"]):
                    continue

            # 3. Instantiate model
            try:
                record = model_cls.model_validate(entity_data)
                loaded_entities[table_name].append(record)
            except Exception as exc:
                raise IngestError(
                    f"Row {row_idx + 1} of {table_config.source_file} failed validation for "
                    f"table {table_name} ({model_cls.__name__}): {exc}\nData: {entity_data}"
                ) from exc

    study = PharmacologyStudy(
        metadata=config.study_metadata,
        subjects=loaded_entities["subjects"],  # type: ignore[arg-type]
        compounds=loaded_entities["compounds"],  # type: ignore[arg-type]
        analytes=loaded_entities["analytes"],  # type: ignore[arg-type]
        assays=loaded_entities["assays"],  # type: ignore[arg-type]
        dose_events=loaded_entities["dosing_events"],  # type: ignore[arg-type]
        biosamples=loaded_entities["biosamples"],  # type: ignore[arg-type]
        observations=loaded_entities["observations"],  # type: ignore[arg-type]
        identity_links=loaded_entities["identity_map"],  # type: ignore[arg-type]
    )
    study.attach_source_files(source_files)

    if validate_schema:
        report = study.check(mode=mode)
        if not report.ok:
            raise ValidationFailed(report)

    return study
