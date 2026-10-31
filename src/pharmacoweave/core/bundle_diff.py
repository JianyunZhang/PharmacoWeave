"""ResultBundle comparative diffing tool.

Compares two ResultBundles:
- Manifest metadata differences (run_id, version, analysis_executed, timestamps)
- File tree changes (added, removed, modified, identical)
- Column and row-level numerical parameter diffs for Parquet tables
- High-level metric diffs for JSON analysis reports (prediction audit, analysis parameters)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import pyarrow.parquet as pq
from pydantic import Field

from pharmacoweave.core.entities import Record


class TableColumnDiff(Record):
    """Numerical or structural differences for one column across two bundles."""

    column_name: str
    diff_type: Literal["numeric", "string", "type_mismatch"]
    max_relative_diff: float | None = None
    max_absolute_diff: float | None = None
    differences_count: int = 0
    sample_diff: str | None = None


class TableDiff(Record):
    """Comparative summary for one tabular artifact across two bundles."""

    relative_path: str
    rows_a: int
    rows_b: int
    columns_only_in_a: list[str] = Field(default_factory=list)
    columns_only_in_b: list[str] = Field(default_factory=list)
    columns_diffs: list[TableColumnDiff] = Field(default_factory=list)


class BundleDiffResult(Record):
    """Structured result of comparing two ResultBundles."""

    bundle_a_path: str
    bundle_b_path: str
    identical: bool
    metadata_diffs: dict[str, tuple[Any, Any]] = Field(default_factory=dict)
    files_added: list[str] = Field(default_factory=list)
    files_removed: list[str] = Field(default_factory=list)
    files_modified: list[str] = Field(default_factory=list)
    files_identical: list[str] = Field(default_factory=list)
    table_diffs: list[TableDiff] = Field(default_factory=list)
    json_diffs: dict[str, Any] = Field(default_factory=dict)
    research_use_only: Literal[True] = True

    def summary_lines(self) -> list[str]:
        """Format human-readable CLI summary lines."""
        lines = [
            f"ResultBundle Diff: {self.bundle_a_path} vs {self.bundle_b_path}",
            f"  Result: {'IDENTICAL' if self.identical else 'DIFFERENT'}",
        ]
        if self.metadata_diffs:
            lines.append("  Metadata Differences:")
            for key, (val_a, val_b) in self.metadata_diffs.items():
                lines.append(f"    - {key}: {val_a!r} -> {val_b!r}")

        lines.append(
            f"  Files: {len(self.files_identical)} identical, "
            f"{len(self.files_modified)} modified, "
            f"{len(self.files_added)} added, "
            f"{len(self.files_removed)} removed"
        )
        if self.files_modified:
            lines.append("  Modified Files:")
            for p in self.files_modified:
                lines.append(f"    * {p}")

        if self.table_diffs:
            lines.append("  Tabular Parameter Diffs:")
            for td in self.table_diffs:
                lines.append(f"    Table: {td.relative_path} (rows: {td.rows_a} vs {td.rows_b})")
                if td.columns_only_in_a:
                    lines.append(f"      - Columns removed in B: {td.columns_only_in_a}")
                if td.columns_only_in_b:
                    lines.append(f"      + Columns added in B: {td.columns_only_in_b}")
                for cd in td.columns_diffs:
                    if cd.max_relative_diff is not None:
                        lines.append(
                            f"      ~ {cd.column_name}: max rel diff = {cd.max_relative_diff:.4e} "
                            f"(abs = {cd.max_absolute_diff:.4e}, {cd.differences_count} diffs)"
                        )
                    else:
                        lines.append(
                            f"      ~ {cd.column_name}: {cd.differences_count} diffs ({cd.sample_diff})"
                        )

        if self.json_diffs:
            lines.append("  JSON Report Differences:")
            for file_key, d in self.json_diffs.items():
                lines.append(f"    File: {file_key}")
                for k, v in d.items():
                    lines.append(f"      - {k}: {v}")
        return lines


def _diff_parquet(
    path_a: Path,
    path_b: Path,
    relative_path: str,
    rel_tol: float,
) -> TableDiff:
    """Compare two Parquet tables and report numerical differences."""
    table_a = pq.read_table(path_a)
    table_b = pq.read_table(path_b)

    cols_a = set(table_a.column_names)
    cols_b = set(table_b.column_names)
    common_cols = sorted(cols_a & cols_b)

    rows_a = table_a.num_rows
    rows_b = table_b.num_rows

    columns_diffs: list[TableColumnDiff] = []
    min_rows = min(rows_a, rows_b)

    for col in common_cols:
        series_a = table_a[col].to_pylist()[:min_rows]
        series_b = table_b[col].to_pylist()[:min_rows]

        # Check types
        sample_val = next((v for v in series_a if v is not None), None)
        if sample_val is None:
            sample_val = next((v for v in series_b if v is not None), None)

        is_numeric = isinstance(sample_val, (int, float)) and not isinstance(sample_val, bool)

        diffs_count = 0
        max_abs = 0.0
        max_rel = 0.0
        sample_diff_str: str | None = None

        for val_a, val_b in zip(series_a, series_b, strict=False):
            if val_a == val_b:
                continue
            if val_a is None or val_b is None:
                diffs_count += 1
                if sample_diff_str is None:
                    sample_diff_str = f"{val_a} vs {val_b}"
                continue

            if is_numeric:
                try:
                    f_a = float(val_a)
                    f_b = float(val_b)
                    abs_diff = abs(f_b - f_a)
                    denominator = max(abs(f_a), abs(f_b), 1e-12)
                    rel_diff = abs_diff / denominator
                    if rel_diff > rel_tol:
                        diffs_count += 1
                        if abs_diff > max_abs:
                            max_abs = abs_diff
                        if rel_diff > max_rel:
                            max_rel = rel_diff
                            sample_diff_str = f"{f_a} vs {f_b}"
                except (ValueError, TypeError):
                    diffs_count += 1
            else:
                diffs_count += 1
                if sample_diff_str is None:
                    sample_diff_str = f"{val_a!r} vs {val_b!r}"

        if diffs_count > 0 or rows_a != rows_b:
            if is_numeric:
                columns_diffs.append(
                    TableColumnDiff(
                        column_name=col,
                        diff_type="numeric",
                        max_relative_diff=max_rel if diffs_count > 0 else 0.0,
                        max_absolute_diff=max_abs if diffs_count > 0 else 0.0,
                        differences_count=diffs_count,
                        sample_diff=sample_diff_str,
                    )
                )
            else:
                columns_diffs.append(
                    TableColumnDiff(
                        column_name=col,
                        diff_type="string",
                        differences_count=diffs_count,
                        sample_diff=sample_diff_str,
                    )
                )

    return TableDiff(
        relative_path=relative_path,
        rows_a=rows_a,
        rows_b=rows_b,
        columns_only_in_a=sorted(cols_a - cols_b),
        columns_only_in_b=sorted(cols_b - cols_a),
        columns_diffs=columns_diffs,
    )


def diff_bundles(
    bundle_a: str | Path,
    bundle_b: str | Path,
    *,
    rel_tol: float = 1e-5,
) -> BundleDiffResult:
    """Compare two ResultBundle directories and produce a structured diff."""
    dir_a = Path(bundle_a)
    dir_b = Path(bundle_b)

    manifest_a_path = dir_a / "manifest.json"
    manifest_b_path = dir_b / "manifest.json"

    if not manifest_a_path.is_file():
        raise FileNotFoundError(f"manifest.json missing in bundle A: {dir_a}")
    if not manifest_b_path.is_file():
        raise FileNotFoundError(f"manifest.json missing in bundle B: {dir_b}")

    raw_a = json.loads(manifest_a_path.read_text(encoding="utf-8"))
    raw_b = json.loads(manifest_b_path.read_text(encoding="utf-8"))

    # Compare metadata
    metadata_diffs: dict[str, tuple[Any, Any]] = {}
    for meta_key in ("study_id", "run_id", "schema_version", "analysis_executed"):
        val_a = raw_a.get(meta_key)
        val_b = raw_b.get(meta_key)
        if val_a != val_b:
            metadata_diffs[meta_key] = (val_a, val_b)

    prod_a = raw_a.get("producer", {}).get("version")
    prod_b = raw_b.get("producer", {}).get("version")
    if prod_a != prod_b:
        metadata_diffs["producer_version"] = (prod_a, prod_b)

    # Compare files list
    files_a = {f["path"].replace("\\", "/"): f["sha256"] for f in raw_a.get("files", [])}
    files_b = {f["path"].replace("\\", "/"): f["sha256"] for f in raw_b.get("files", [])}

    disk_files_a = {
        str(f.relative_to(dir_a)).replace("\\", "/")
        for f in dir_a.rglob("*")
        if f.is_file() and f.name not in {"manifest.sha256"}
    }
    disk_files_b = {
        str(f.relative_to(dir_b)).replace("\\", "/")
        for f in dir_b.rglob("*")
        if f.is_file() and f.name not in {"manifest.sha256"}
    }

    all_keys_a = set(files_a.keys()) | disk_files_a
    all_keys_b = set(files_b.keys()) | disk_files_b

    added = sorted(all_keys_b - all_keys_a)
    removed = sorted(all_keys_a - all_keys_b)
    common = sorted(all_keys_a & all_keys_b)

    identical: list[str] = []
    modified: list[str] = []

    for p in common:
        disk_a = dir_a / p
        disk_b = dir_b / p
        sha_a = files_a.get(p)
        sha_b = files_b.get(p)
        if disk_a.is_file() and disk_b.is_file():
            # Quick size check before reading
            if sha_a is not None and sha_b is not None and sha_a == sha_b:
                if disk_a.stat().st_size == disk_b.stat().st_size:
                    identical.append(p)
                else:
                    modified.append(p)
            elif disk_a.read_bytes() == disk_b.read_bytes():
                identical.append(p)
            else:
                modified.append(p)
        elif disk_a.is_file() or disk_b.is_file():
            # Exists on disk in one bundle but missing in the other
            modified.append(p)
        elif sha_a == sha_b:
            identical.append(p)
        else:
            modified.append(p)

    # Compare tables for modified parquet files
    table_diffs: list[TableDiff] = []
    json_diffs: dict[str, Any] = {}

    for mod_path in modified:
        disk_a = dir_a / mod_path
        disk_b = dir_b / mod_path
        if not disk_a.is_file() or not disk_b.is_file():
            continue

        if mod_path.endswith(".parquet"):
            try:
                td = _diff_parquet(disk_a, disk_b, mod_path, rel_tol=rel_tol)
                table_diffs.append(td)
            except Exception as exc:
                table_diffs.append(
                    TableDiff(
                        relative_path=mod_path,
                        rows_a=-1,
                        rows_b=-1,
                        columns_diffs=[
                            TableColumnDiff(
                                column_name="_error",
                                diff_type="string",
                                sample_diff=f"Failed to read parquet: {exc}",
                            )
                        ],
                    )
                )
        elif mod_path.endswith(".json") and mod_path != "manifest.json":
            try:
                obj_a = json.loads(disk_a.read_text(encoding="utf-8"))
                obj_b = json.loads(disk_b.read_text(encoding="utf-8"))
                if isinstance(obj_a, dict) and isinstance(obj_b, dict):
                    diff_keys = {
                        k
                        for k in set(obj_a.keys()) | set(obj_b.keys())
                        if obj_a.get(k) != obj_b.get(k)
                    }
                    if diff_keys:
                        json_diffs[mod_path] = {
                            k: f"{obj_a.get(k)!r} -> {obj_b.get(k)!r}" for k in sorted(diff_keys)
                        }
            except Exception:
                pass

    is_identical = not (metadata_diffs or added or removed or modified)

    return BundleDiffResult(
        bundle_a_path=str(dir_a),
        bundle_b_path=str(dir_b),
        identical=is_identical,
        metadata_diffs=metadata_diffs,
        files_added=added,
        files_removed=removed,
        files_modified=modified,
        files_identical=identical,
        table_diffs=table_diffs,
        json_diffs=json_diffs,
    )
