"""Result bundles and the command line."""

import hashlib
import json
from pathlib import Path

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
    oral_cfg = load_nca_config(repo_root / "configs" / "nca_single_dose_oral.yaml")
    assert oral_cfg.design_route == "oral"
    inf_cfg = load_nca_config(repo_root / "configs" / "nca_single_dose_iv_infusion.yaml")
    assert inf_cfg.design_route == "iv_infusion"

    incomplete = tmp_path / "incomplete.yaml"
    incomplete.write_text("schema_version: '0.1.0'\nanalysis: nca\ndesign_route: oral\n", encoding="utf-8")
    with pytest.raises(IngestError, match="implemented NCA design"):
        load_nca_config(incomplete)

    multi_dose = tmp_path / "multiple_dose.yaml"
    multi_dose.write_text(
        (repo_root / "configs" / "nca_single_dose_oral.yaml")
        .read_text(encoding="utf-8")
        .replace("dosing: single_dose", "dosing: multiple_dose"),
        encoding="utf-8",
    )
    with pytest.raises(IngestError, match="implemented NCA design"):
        load_nca_config(multi_dose)

    sparse = tmp_path / "sparse.yaml"
    sparse.write_text(
        (repo_root / "configs" / "nca_single_dose_oral.yaml")
        .read_text(encoding="utf-8")
        .replace("sampling: rich_serial", "sampling: sparse"),
        encoding="utf-8",
    )
    with pytest.raises(IngestError, match="implemented NCA design"):
        load_nca_config(sparse)


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

    oral_dir = repo_root / "examples" / "data" / "oral_single_dose"
    oral_bundle = tmp_path / "cli-oral-bundle"
    assert (
        main(
            [
                "nca",
                str(oral_dir),
                "--config",
                str(repo_root / "configs" / "nca_single_dose_oral.yaml"),
                "--with-fit",
                "--bundle",
                str(oral_bundle),
                "--run-id",
                "cli-oral",
            ]
        )
        == 0
    )
    captured_oral = capsys.readouterr().out
    assert "CL/F=" in captured_oral and "Vz/F=" in captured_oral
    assert (oral_bundle / "analysis" / "nca_parameters.parquet").is_file()
    assert (oral_bundle / "analysis" / "model_parameters.parquet").is_file()
    oral_html = (oral_bundle / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "CL/F=" in oral_html and "one_compartment_extravascular" in oral_html

    two_comp_dir = repo_root / "examples" / "data" / "two_compartment_iv_bolus"
    two_comp_bundle = tmp_path / "cli-2c-bundle"
    assert (
        main(
            [
                "fit",
                str(two_comp_dir),
                "--model",
                "two_compartment_iv_bolus",
                "--bundle",
                str(two_comp_bundle),
                "--run-id",
                "cli-2c",
            ]
        )
        == 0
    )
    captured_2c = capsys.readouterr().out
    assert "V1=" in captured_2c and "Q=" in captured_2c and "V2=" in captured_2c
    assert (two_comp_bundle / "analysis" / "model_parameters.parquet").is_file()


def test_bundle_persists_ivive_and_prediction_audit(tmp_path, repo_root) -> None:
    from pharmacoweave.diagnostics.prediction_audit import audit_study_predictions
    from pharmacoweave.dmpk import (
        AdmeAssayKind,
        AdmeStatus,
        InVitroAdmeRecord,
        PhysiologicalScalingParams,
        predict_hepatic_clearance_well_stirred,
    )
    from pharmacoweave.integrations.pharmacology_evidence import collect_run_evidence

    study = build_study()
    nca = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    clint_rec = InVitroAdmeRecord(
        record_id="ADME-MIC-CLI",
        compound_id="PW-0001",
        assay_kind=AdmeAssayKind.MICROSOMAL_CLEARANCE,
        system="rat liver microsomes",
        status=AdmeStatus.MEASURED,
        limitations=["Research use only. Not a dosing recommendation."],
        value=30.0,
        unit="uL/min/mg_protein",
    )
    phys = PhysiologicalScalingParams(
        profile_id="PHYS-RAT",
        species="Rattus norvegicus",
        body_weight_kg=0.25,
        liver_weight_g=10.0,
        hepatic_blood_flow_l_h=1.0,
        microsomal_protein_mg_per_g_liver=45.0,
        source="Synthetic test fixture",
    )
    ivive_pred = predict_hepatic_clearance_well_stirred(
        clint_rec,
        physiology=phys,
        fu_plasma=0.1,
        fu_inc=1.0,
        blood_to_plasma_ratio=1.0,
    )

    estimates = collect_run_evidence(run_id="run-nca", nca_parameters=nca)
    audit_rep = audit_study_predictions(
        study,
        estimates,
        priors=[
            from_ivive
            for from_ivive in collect_run_evidence(run_id="run-ivive", ivive_predictions=[ivive_pred])
            if from_ivive.evidence_type.value == "predicted"
        ],
        audit_id="AUDIT-TEST",
    )

    destination = tmp_path / "ivive-audit-bundle"
    manifest = study.export_bundle(
        destination,
        run_id="ivive-audit-run",
        nca_parameters=nca,
        ivive_predictions=[ivive_pred],
        prediction_audit=audit_rep,
    )
    assert manifest.analysis_executed is True
    assert manifest.data["ivive_predictions"] == "analysis/ivive_predictions.parquet"
    assert manifest.data["prediction_audit"] == "analysis/prediction_audit.json"

    ivive_table = pq.read_table(destination / "analysis" / "ivive_predictions.parquet")
    assert ivive_table.num_rows == 1
    assert ivive_table.column("species").to_pylist() == ["Rattus norvegicus"]
    assert ivive_table.column("evidence_grade").to_pylist() == ["predicted"]

    audit_payload = json.loads((destination / "analysis" / "prediction_audit.json").read_text(encoding="utf-8"))
    assert audit_payload["audit_id"] == "AUDIT-TEST"
    assert len(audit_payload["entries"]) == 2

    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "Exploratory IVIVE (Well-Stirred)" in html
    assert "Prior vs. Posterior Prediction Audit" in html
    assert "RAT" in html

    exchange = EvidenceExchange.model_validate_json(
        (destination / "evidence" / "pharmacology_evidence.json").read_text(encoding="utf-8")
    )
    assert any(
        rec.claim_type == "ivive_hepatic_plasma_clearance" and rec.evidence_type.value == "predicted"
        for rec in exchange.records
    )


def test_cli_ivive_and_audit_predictions(tmp_path, repo_root, capsys) -> None:
    ivive_json = tmp_path / "ivive_input.json"
    ivive_json.write_text(
        json.dumps(
            {
                "clint_record": {
                    "record_id": "ADME-MIC-CLI",
                    "compound_id": "PW-0001",
                    "assay_kind": "microsomal_clearance",
                    "system": "rat liver microsomes",
                    "status": "measured",
                    "limitations": ["Research use only. Not a dosing recommendation."],
                    "value": 30.0,
                    "unit": "uL/min/mg_protein",
                },
                "physiology": {
                    "profile_id": "PHYS-RAT-CLI",
                    "species": "Rattus norvegicus",
                    "body_weight_kg": 0.25,
                    "liver_weight_g": 10.0,
                    "hepatic_blood_flow_l_h": 1.0,
                    "microsomal_protein_mg_per_g_liver": 45.0,
                    "source": "CLI test fixture",
                },
                "fu_plasma": 0.1,
                "fu_inc": 1.0,
                "blood_to_plasma_ratio": 1.0,
            }
        ),
        encoding="utf-8",
    )
    evidence_out = tmp_path / "ivive_evidence.json"
    ivive_bundle = tmp_path / "ivive_cli_bundle"

    assert (
        main(
            [
                "ivive",
                str(ivive_json),
                "--evidence-out",
                str(evidence_out),
                "--study",
                str(EXAMPLE_DIR),
                "--bundle",
                str(ivive_bundle),
                "--run-id",
                "cli-ivive",
            ]
        )
        == 0
    )
    captured = capsys.readouterr().out
    assert "CLH_plasma=" in captured
    assert evidence_out.is_file()
    assert (ivive_bundle / "analysis" / "ivive_predictions.parquet").is_file()

    # CLI audit-predictions
    audit_out = tmp_path / "audit_report.json"
    audit_bundle = tmp_path / "audit_cli_bundle"
    assert (
        main(
            [
                "audit-predictions",
                str(EXAMPLE_DIR),
                "--priors",
                str(evidence_out),
                "--nca-config",
                str(repo_root / "configs" / "nca_single_dose_iv.yaml"),
                "--output",
                str(audit_out),
                "--bundle",
                str(audit_bundle),
                "--run-id",
                "cli-audit",
            ]
        )
        == 0
    )
    captured_audit = capsys.readouterr().out
    assert "clearance" in captured_audit
    assert audit_out.is_file()
    assert (audit_bundle / "analysis" / "prediction_audit.json").is_file()


def test_cli_verify_bundle_and_diff(tmp_path: Path, repo_root: Path, capsys) -> None:
    import shutil

    study = build_study()
    bundle_a = tmp_path / "cli_bundle_a"
    bundle_b = tmp_path / "cli_bundle_b"

    study.export_bundle(bundle_a, run_id="run_a")
    shutil.copytree(bundle_a, bundle_b)

    # 1. verify-bundle command
    assert main(["verify-bundle", str(bundle_a)]) == 0
    captured = capsys.readouterr().out
    assert "Status: PASSED" in captured

    # 2. diff command identical
    assert main(["diff", str(bundle_a), str(bundle_b)]) == 0
    captured_diff = capsys.readouterr().out
    assert "Result: IDENTICAL" in captured_diff

    # 3. verify-bundle with JSON flag
    assert main(["verify-bundle", str(bundle_a), "--json"]) == 0
    json_out = capsys.readouterr().out
    assert '"is_valid": true' in json_out

    # 4. Tamper bundle_b: verify fails and diff detects differences
    (bundle_b / "study" / "compounds.parquet").unlink()
    assert main(["verify-bundle", str(bundle_b)]) == 1
    captured_fail = capsys.readouterr().out
    assert "Status: FAILED" in captured_fail

    assert main(["diff", str(bundle_a), str(bundle_b)]) == 1
    captured_diff_fail = capsys.readouterr().out
    assert "Result: DIFFERENT" in captured_diff_fail


def test_cli_ingest_mapping(tmp_path: Path, capsys) -> None:
    from tests.unit.test_ingest_mapping import _setup_raw_data

    mapping_path = _setup_raw_data(tmp_path)
    output_dir = tmp_path / "ingested_study"
    bundle_dir = tmp_path / "ingested_bundle"

    assert (
        main(
            [
                "ingest",
                "--mapping",
                str(mapping_path),
                "--output-dir",
                str(output_dir),
                "--bundle",
                str(bundle_dir),
                "--run-id",
                "run-ingest-cli",
            ]
        )
        == 0
    )
    captured = capsys.readouterr().out
    assert "Wrote canonical study directory" in captured
    assert "Wrote result bundle" in captured

    assert (output_dir / "study_metadata.json").is_file()
    assert (output_dir / "subjects.csv").is_file()
    assert (bundle_dir / "manifest.json").is_file()


