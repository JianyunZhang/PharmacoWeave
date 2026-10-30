"""Result bundles and the command line."""

import hashlib
import json

import pyarrow.parquet as pq
import pytest
from examples.direct_effect_and_evidence import E0, build_adme, build_effects, build_exposures, build_tox
from examples.iv_bolus_case import EXAMPLE_DIR, build_study

from pharmacoweave.cli import main
from pharmacoweave.core.evidence import EvidenceExchange
from pharmacoweave.errors import IngestError
from pharmacoweave.pk_nca import load_nca_config
from pharmacoweave.pkpd import fit_linked_direct_emax


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


def test_nca_bundle_stores_parameters_semilog_derived_and_evidence(tmp_path, repo_root) -> None:
    study = build_study()
    results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    destination = tmp_path / "nca-bundle"
    manifest = study.export_bundle(destination, run_id="nca-run", nca_parameters=results)
    assert manifest.analysis_executed is True
    assert manifest.data["nca_parameters"] == "analysis/nca_parameters.parquet"
    assert manifest.data["derived_observations"] == "analysis/derived_observations.parquet"
    assert not any(item.role == "nca_parameters" for item in manifest.omitted)
    table = pq.read_table(destination / "analysis" / "nca_parameters.parquet")
    assert table.num_rows == 2
    derived = pq.read_table(destination / "analysis" / "derived_observations.parquet")
    assert derived.num_rows == len(study.observations)
    assert (destination / "figures" / "concentration_time.svg").is_file()
    assert (destination / "figures" / "concentration_time_semilog.svg").is_file()
    exchange = EvidenceExchange.model_validate_json(
        (destination / "evidence" / "pharmacology_evidence.json").read_text(encoding="utf-8")
    )
    assert any(rec.claim_type == "auc_inf" and rec.evidence_type.value == "estimated" for rec in exchange.records)
    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "NCA was run" in html
    assert "were not run" not in html


def test_phase2_bundle_persists_emax_adme_tox_and_screened_evidence(tmp_path, repo_root) -> None:
    study = build_study()
    nca = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    _linked, emax_fit = fit_linked_direct_emax(build_effects(), build_exposures(), e0=E0, e0_unit="percent")
    adme = build_adme()
    tox = [build_tox()]

    destination = tmp_path / "phase2-bundle"
    manifest = study.export_bundle(
        destination,
        run_id="phase2-run",
        nca_parameters=nca,
        emax_fits=[emax_fit],
        adme_records=adme,
        tox_endpoints=tox,
    )
    assert manifest.analysis_executed is True
    assert manifest.data["emax_parameters"] == "analysis/emax_parameters.parquet"
    assert manifest.data["emax_predictions"] == "analysis/emax_predictions.parquet"
    assert manifest.data["in_vitro_adme"] == "experimental/in_vitro_adme.parquet"
    assert manifest.data["tox_endpoints"] == "experimental/tox_endpoints.parquet"
    assert manifest.data["concentration_effect_figure"] == "figures/concentration_effect.svg"

    emax_params = pq.read_table(destination / "analysis" / "emax_parameters.parquet")
    assert set(emax_params.column("name").to_pylist()) == {"E0", "Emax", "EC50", "hill"}
    emax_preds = pq.read_table(destination / "analysis" / "emax_predictions.parquet")
    assert set(emax_preds.column("evidence_grade").to_pylist()) == {"predicted"}
    adme_table = pq.read_table(destination / "experimental" / "in_vitro_adme.parquet")
    assert adme_table.num_rows == len(adme)
    tox_table = pq.read_table(destination / "experimental" / "tox_endpoints.parquet")
    assert tox_table.num_rows == 1

    exchange = EvidenceExchange.model_validate_json(
        (destination / "evidence" / "pharmacology_evidence.json").read_text(encoding="utf-8")
    )
    grades = {(rec.entity_type, rec.evidence_type.value) for rec in exchange.records}
    assert ("NcaParameterSet", "estimated") in grades
    assert ("DirectEmaxFit", "estimated") in grades
    assert ("PredictedEffect", "predicted") in grades
    assert ("InVitroAdmeRecord", "observed") in grades
    assert ("ToxEndpoint", "observed") in grades


def test_nca_config_rejects_designs_outside_the_locked_scope(tmp_path, repo_root) -> None:
    loaded = load_nca_config(repo_root / "configs" / "nca_single_dose_iv.yaml")
    assert loaded.design_route == "iv_bolus"
    assert loaded.blq_rule == "keep_as_censored"
    path = tmp_path / "oral.yaml"
    path.write_text("schema_version: '0.1.0'\nanalysis: nca\ndesign_route: oral\n", encoding="utf-8")
    with pytest.raises(IngestError, match="implemented NCA design"):
        load_nca_config(path)


def test_cli_validates_decouples_nca_fits_wnls_and_runs_emax(tmp_path, repo_root, capsys) -> None:
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
    for name in ("pd_models", "pkpd", "dmpk", "toxicology"):
        line = next(item for item in text.splitlines() if item.startswith(name))
        assert "partial" in line
        assert "reserved" not in line
    assert main(["validate", str(repo_root / "examples")]) == 2

    nca_bundle = tmp_path / "cli-nca-pure"
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
                "cli-nca-pure",
            ]
        )
        == 0
    )
    assert (nca_bundle / "analysis" / "nca_parameters.parquet").is_file()
    assert (nca_bundle / "figures" / "concentration_time.svg").is_file()
    assert (nca_bundle / "figures" / "concentration_time_semilog.svg").is_file()
    assert not (nca_bundle / "figures" / "goodness_of_fit.svg").exists()

    nca_with_fit = tmp_path / "cli-nca-with-fit"
    assert (
        main(
            [
                "nca",
                str(EXAMPLE_DIR),
                "--config",
                str(repo_root / "configs" / "nca_single_dose_iv.yaml"),
                "--with-fit",
                "--bundle",
                str(nca_with_fit),
                "--run-id",
                "cli-nca-with-fit",
            ]
        )
        == 0
    )
    assert (nca_with_fit / "figures" / "goodness_of_fit.svg").is_file()
    gof = (nca_with_fit / "figures" / "goodness_of_fit.svg").read_text(encoding="utf-8")
    assert gof.count("<circle") == 16
    assert "not plotted as zero" in gof
    capsys.readouterr()

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
    assert (fit_bundle / "figures" / "weighted_residuals.svg").is_file()
    captured_fit = capsys.readouterr().out
    assert "CL=" in captured_fit
    assert "status=fail" not in captured_fit

    emax_input = repo_root / "examples" / "data" / "direct_effect_emax.json"
    evidence_out = tmp_path / "emax_evidence.json"
    emax_bundle = tmp_path / "cli-emax-bundle"
    assert (
        main(
            [
                "emax",
                str(emax_input),
                "--evidence-out",
                str(evidence_out),
                "--study",
                str(EXAMPLE_DIR),
                "--bundle",
                str(emax_bundle),
                "--run-id",
                "cli-emax",
            ]
        )
        == 0
    )
    captured_emax = capsys.readouterr().out
    assert "method=direct_emax_nls_v1" in captured_emax
    assert "EC50=" in captured_emax
    assert evidence_out.is_file()
    assert (emax_bundle / "analysis" / "emax_parameters.parquet").is_file()
    assert (emax_bundle / "figures" / "concentration_effect.svg").is_file()
