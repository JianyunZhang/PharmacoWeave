"""Result bundles and the command line."""

import hashlib
import json

import pyarrow.parquet as pq
import pytest
from examples.iv_bolus_case import EXAMPLE_DIR, build_study

from pharmacoweave.cli import main
from pharmacoweave.errors import IngestError
from pharmacoweave.pk_nca import load_nca_config


def test_bundle_is_research_use_and_does_not_invent_nca(tmp_path, repo_root) -> None:
    study = build_study()
    destination = tmp_path / "bundle"
    manifest = study.export_bundle(destination, run_id="run-001")
    assert manifest.quality_status == "research_only"
    assert manifest.analysis_executed is False
    assert "nca_parameters" not in manifest.data
    assert any(item.role == "nca_parameters" for item in manifest.omitted)
    manifest_bytes = (destination / "manifest.json").read_bytes()
    digest, _name = (destination / "manifest.sha256").read_text(encoding="utf-8").split()
    assert digest == hashlib.sha256(manifest_bytes).hexdigest()
    for item in manifest.files:
        found = hashlib.sha256((destination / item.path).read_bytes()).hexdigest()
        assert found == item.sha256
    concentrations = pq.read_table(destination / "experimental" / "concentration_observations.parquet")
    effects = pq.read_table(destination / "experimental" / "pd_observations.parquet")
    assert concentrations.num_rows == 20
    assert effects.num_rows == 0
    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "Research use only" in html
    assert "were not run" in html
    assert "PW-SYNTH-IV-001" in html
    payload = json.loads(manifest_bytes)
    assert payload["schema_name"] == "pharmacology-result-bundle"
    with pytest.raises(IngestError):
        study.export_bundle(destination, run_id="run-002")
    nca = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    assert len(nca.nca_parameters) == 2
    assert all(item.auc_last is not None for item in nca.nca_parameters)


def test_nca_config_rejects_designs_outside_the_locked_scope(tmp_path, repo_root) -> None:
    loaded = load_nca_config(repo_root / "configs" / "nca_single_dose_iv.yaml")
    assert loaded.design_route == "iv_bolus"
    assert loaded.blq_rule == "keep_as_censored"
    path = tmp_path / "oral.yaml"
    path.write_text("schema_version: '0.1.0'\nanalysis: nca\ndesign_route: oral\n", encoding="utf-8")
    with pytest.raises(IngestError, match="implemented NCA design"):
        load_nca_config(path)


def test_cli_validates_the_example_and_reports_capabilities(tmp_path, repo_root, capsys) -> None:
    assert main(["validate", str(EXAMPLE_DIR), "--strict"]) == 0
    captured = capsys.readouterr()
    assert "passed" in captured.out
    bundle = tmp_path / "cli-bundle"
    assert main(["validate", str(EXAMPLE_DIR), "--bundle", str(bundle), "--run-id", "cli-run"]) == 0
    assert (bundle / "manifest.json").is_file()
    assert main(["capabilities"]) == 0
    text = capsys.readouterr().out
    rows = {line.split()[0]: line for line in text.splitlines() if line.strip()}
    assert "available" in rows["pk_nca"]
    assert "available" in rows["pk_models"]
    assert "partial" in rows["integrations"]
    assert "reserved" in rows["pbpk"]
    assert main(["validate", str(repo_root / "examples")]) == 2
