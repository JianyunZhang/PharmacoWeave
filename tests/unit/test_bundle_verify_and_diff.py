"""Unit tests for offline ResultBundle verification (verify-bundle) and diffing (diff)."""

import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.bundle_diff import diff_bundles
from pharmacoweave.core.bundle_verify import verify_bundle


def test_verify_bundle_clean_passes(tmp_path: Path, repo_root: Path) -> None:
    study = build_study()
    nca_results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    bundle_dir = tmp_path / "clean_bundle"
    study.export_bundle(bundle_dir, run_id="run_clean", nca_parameters=nca_results)

    result = verify_bundle(bundle_dir)
    assert result.is_valid is True
    assert result.run_id == "run_clean"
    assert result.study_id == "PW-SYNTH-IV-001"
    assert result.files_checked > 10
    assert result.files_ok == result.files_checked
    assert len(result.errors) == 0

    lines = result.summary_lines()
    assert any("Status: PASSED" in line for line in lines)


def test_verify_bundle_tampered_fails(tmp_path: Path) -> None:
    study = build_study()
    bundle_dir = tmp_path / "tampered_bundle"
    study.export_bundle(bundle_dir, run_id="run_tampered")

    # Tamper with a file by appending a byte
    target_file = bundle_dir / "study" / "compounds.parquet"
    target_file.write_bytes(target_file.read_bytes() + b"\x00")

    result = verify_bundle(bundle_dir)
    assert result.is_valid is False
    assert any("Checksum mismatch" in err and "compounds.parquet" in err for err in result.errors)


def test_verify_bundle_missing_file_fails(tmp_path: Path) -> None:
    study = build_study()
    bundle_dir = tmp_path / "missing_file_bundle"
    study.export_bundle(bundle_dir, run_id="run_missing")

    # Remove subjects.parquet
    (bundle_dir / "study" / "subjects.parquet").unlink()

    result = verify_bundle(bundle_dir)
    assert result.is_valid is False
    assert any("Missing file declared in manifest" in err for err in result.errors)


def test_diff_bundles_identical(tmp_path: Path) -> None:
    study = build_study()
    bundle_a = tmp_path / "bundle_a"
    bundle_b = tmp_path / "bundle_b"

    study.export_bundle(bundle_a, run_id="run_test")
    shutil.copytree(bundle_a, bundle_b)

    diff_result = diff_bundles(bundle_a, bundle_b)
    assert diff_result.identical is True
    assert len(diff_result.files_modified) == 0
    assert len(diff_result.files_added) == 0
    assert len(diff_result.files_removed) == 0
    assert len(diff_result.metadata_diffs) == 0

    lines = diff_result.summary_lines()
    assert any("Result: IDENTICAL" in line for line in lines)


def test_diff_bundles_with_changes(tmp_path: Path, repo_root: Path) -> None:
    study = build_study()
    nca_results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")

    bundle_a = tmp_path / "bundle_orig"
    bundle_b = tmp_path / "bundle_modified"

    study.export_bundle(bundle_a, run_id="run_v1", nca_parameters=nca_results)
    shutil.copytree(bundle_a, bundle_b)

    # 1. Change run_id in bundle_b manifest
    manifest_b_file = bundle_b / "manifest.json"
    manifest_b_data = json.loads(manifest_b_file.read_text(encoding="utf-8"))
    manifest_b_data["run_id"] = "run_v2"
    manifest_b_file.write_text(json.dumps(manifest_b_data, indent=2), encoding="utf-8")

    # 2. Modify numeric value in nca_parameters.parquet of bundle_b
    nca_table_path = bundle_b / "analysis" / "nca_parameters.parquet"
    table = pq.read_table(nca_table_path)
    # Alter cmax column
    pydict = table.to_pydict()
    pydict["cmax"][0] = pydict["cmax"][0] * 1.50  # 50% increase
    new_table = pa.Table.from_pydict(pydict)
    pq.write_table(new_table, nca_table_path)

    # 3. Add an extra file in bundle_b
    (bundle_b / "extra_note.txt").write_text("extra", encoding="utf-8")

    diff_result = diff_bundles(bundle_a, bundle_b)
    assert diff_result.identical is False
    assert "run_id" in diff_result.metadata_diffs
    assert diff_result.metadata_diffs["run_id"] == ("run_v1", "run_v2")

    # Check table diff
    nca_diff = next(
        (td for td in diff_result.table_diffs if "nca_parameters.parquet" in td.relative_path),
        None,
    )
    assert nca_diff is not None
    cmax_diff = next((cd for cd in nca_diff.columns_diffs if cd.column_name == "cmax"), None)
    assert cmax_diff is not None
    assert cmax_diff.max_relative_diff == pytest.approx(1.0 / 3.0, rel=1e-2)  # abs(1.5 - 1.0)/1.5 = 0.3333

    lines = diff_result.summary_lines()
    assert any("Result: DIFFERENT" in line for line in lines)
    assert any("run_id: 'run_v1' -> 'run_v2'" in line for line in lines)
