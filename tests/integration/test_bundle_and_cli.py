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
    assert len(nca) == 2
    assert nca[0].auc_inf == pytest.approx(5.0, rel=1e-4)


def test_nca_bundle_stores_parameters_and_a_curve(tmp_path, repo_root) -> None:
    study = build_study()
    results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    destination = tmp_path / "nca-bundle"
    manifest = study.export_bundle(destination, run_id="nca-run", nca_parameters=results)
    assert manifest.analysis_executed is True
    assert manifest.data["nca_parameters"] == "analysis/nca_parameters.parquet"
    assert not any(item.role == "nca_parameters" for item in manifest.omitted)
    table = pq.read_table(destination / "analysis" / "nca_parameters.parquet")
    assert table.num_rows == 2
    assert (destination / "figures" / "concentration_time.svg").is_file()
    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "NCA was run" in html
    assert "were not run" not in html


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
    assert "pk_nca" in text and "partial" in text
    assert "pbpk" in text and "reserved" in text
    assert main(["validate", str(repo_root / "examples")]) == 2
    nca_bundle = tmp_path / "cli-nca"
    assert (
        main(
            [
                "nca",
                str(EXAMPLE_DIR),
                "--config",
                str(repo_root / "configs" / "nca_single_dose_iv.yaml"),
                "--bundle",
                str(nca_bundle),
                "--run-id",
                "cli-nca",
            ]
        )
        == 0
    )
    assert (nca_bundle / "analysis" / "nca_parameters.parquet").is_file()
    assert (nca_bundle / "figures" / "goodness_of_fit.svg").is_file()
    assert (nca_bundle / "figures" / "residual_diagnostics.svg").is_file()
    gof = (nca_bundle / "figures" / "goodness_of_fit.svg").read_text(encoding="utf-8")
    assert gof.count("<circle") == 16
    assert "not plotted as zero" in gof
    fit_bundle = tmp_path / "cli-fit"
    assert (
        main(
            [
                "fit",
                str(EXAMPLE_DIR),
                "--method",
                "weighted_nls",
                "--weight",
                "1/y^2",
                "--bundle",
                str(fit_bundle),
                "--run-id",
                "cli-fit",
            ]
        )
        == 0
    )
    assert (fit_bundle / "figures" / "residual_diagnostics.svg").is_file()
    captured_fit = capsys.readouterr().out
    assert "CL=" in captured_fit
    assert "status=fail" not in captured_fit
