"""ResultBundle offline cryptographic and schema verification.

Verifies the integrity of a written ResultBundle:
1. Validates manifest.json against the ResultBundle schema.
2. Re-hashes every file listed in the manifest and checks against its recorded SHA-256 digest.
3. Detects missing files, corrupted files, and extraneous unmanifested files.
4. Checks provenance and QC validation records.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Literal

import jsonschema
from pydantic import Field

from pharmacoweave.core.entities import Record
from pharmacoweave.core.provenance import ResultBundleManifest
from pharmacoweave.schema_export import schema_directory


class BundleVerificationResult(Record):
    """Structured report returned by bundle verification."""

    bundle_path: str
    is_valid: bool
    manifest_version: str | None = None
    run_id: str | None = None
    study_id: str | None = None
    files_checked: int = 0
    files_ok: int = 0
    errors: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    checked_files: dict[str, str] = Field(default_factory=dict)
    research_use_only: Literal[True] = True

    def summary_lines(self) -> list[str]:
        """Format a human-readable CLI summary."""
        lines = [
            f"ResultBundle Verification: {self.bundle_path}",
            f"  Status: {'PASSED' if self.is_valid else 'FAILED'}",
            f"  Run ID: {self.run_id or 'unknown'}",
            f"  Study ID: {self.study_id or 'unknown'}",
            f"  Manifest Version: {self.manifest_version or 'unknown'}",
            f"  Files Checked: {self.files_checked} ({self.files_ok} valid)",
        ]
        if self.errors:
            lines.append("  Errors:")
            for err in self.errors:
                lines.append(f"    - {err}")
        if self.warnings:
            lines.append("  Warnings:")
            for warn in self.warnings:
                lines.append(f"    - {warn}")
        return lines


def _sha256(path: Path) -> str:
    """Calculate the hex SHA-256 digest of a file."""
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def verify_bundle(bundle_dir: str | Path) -> BundleVerificationResult:
    """Perform an offline verification of a ResultBundle directory.

    Checks:
    - manifest.json presence and schema compliance
    - SHA-256 checksum of every recorded file
    - Provenance and QC records
    - Absence of unmanifested artifacts
    """
    directory = Path(bundle_dir)
    errors: list[str] = []
    warnings: list[str] = []
    checked_files: dict[str, str] = {}

    if not directory.is_dir():
        return BundleVerificationResult(
            bundle_path=str(directory),
            is_valid=False,
            errors=[f"Target path {directory} is not a directory."],
        )

    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        return BundleVerificationResult(
            bundle_path=str(directory),
            is_valid=False,
            errors=[f"manifest.json is missing in {directory}."],
        )

    try:
        manifest_raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        return BundleVerificationResult(
            bundle_path=str(directory),
            is_valid=False,
            errors=[f"manifest.json is not valid JSON: {exc}"],
        )

    # 1. Schema validation against checked-in JSON schema
    try:
        schema_path = schema_directory() / "result_bundle.schema.json"
        if schema_path.is_file():
            bundle_schema = json.loads(schema_path.read_text(encoding="utf-8"))
            jsonschema.validate(instance=manifest_raw, schema=bundle_schema)
    except Exception as exc:
        errors.append(f"manifest.json failed result_bundle schema validation: {exc}")

    # 2. Pydantic model validation
    try:
        manifest = ResultBundleManifest.model_validate(manifest_raw)
    except Exception as exc:
        errors.append(f"manifest.json failed ResultBundleManifest validation: {exc}")
        return BundleVerificationResult(
            bundle_path=str(directory),
            is_valid=False,
            errors=errors,
            warnings=warnings,
        )

    if manifest.quality_status != "research_only":
        errors.append(f"quality_status must be 'research_only', got {manifest.quality_status!r}.")

    # 3. Check file integrity
    manifested_relative_paths: set[str] = {"manifest.json"}
    files_ok = 0

    for item in manifest.files:
        manifested_relative_paths.add(item.path.replace("\\", "/"))
        file_disk_path = directory / item.path
        if not file_disk_path.is_file():
            errors.append(f"Missing file declared in manifest: {item.path}")
            checked_files[item.path] = "missing"
            continue

        try:
            actual_sha = _sha256(file_disk_path)
            if actual_sha != item.sha256:
                errors.append(
                    f"Checksum mismatch for {item.path}: expected {item.sha256}, got {actual_sha}."
                )
                checked_files[item.path] = "hash_mismatch"
            else:
                files_ok += 1
                checked_files[item.path] = "ok"
        except Exception as exc:
            errors.append(f"Error reading file {item.path}: {exc}")
            checked_files[item.path] = "read_error"

    # 4. Check for extraneous unmanifested files
    for found_path in directory.rglob("*"):
        if found_path.is_file():
            rel = str(found_path.relative_to(directory)).replace("\\", "/")
            if rel not in manifested_relative_paths:
                warnings.append(f"Unmanifested extraneous file present: {rel}")

    # 5. Check provenance source files
    prov_path = directory / manifest.provenance
    if not prov_path.is_file():
        errors.append(f"Provenance file {manifest.provenance} missing.")

    # 6. Check QC report
    qc_path = directory / "qc" / "ingest_validation.json"
    if not qc_path.is_file():
        warnings.append("qc/ingest_validation.json missing.")

    is_valid = len(errors) == 0

    return BundleVerificationResult(
        bundle_path=str(directory),
        is_valid=is_valid,
        manifest_version=manifest.schema_version,
        run_id=manifest.run_id,
        study_id=manifest.study_id,
        files_checked=len(manifest.files),
        files_ok=files_ok,
        errors=errors,
        warnings=warnings,
        checked_files=checked_files,
    )
