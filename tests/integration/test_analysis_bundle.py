"""An executed analysis writes figures and sets analysis_executed."""

import json

import pyarrow.parquet as pq
from examples.iv_bolus_case import EXAMPLE_DIR, build_study

from pharmacoweave.cli import main
from pharmacoweave.core.analysis_export import AnalysisPayload
from pharmacoweave.core.bundle import export_bundle
from pharmacoweave.core.entities import Route
from pharmacoweave.core.modeling import (
    AnalysisRecords,
    DiagnosticRecord,
    DiagnosticStatus,
    ModelFit,
    NcaParameterSet,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.version import __version__


def test_analysis_bundle_records_nca_fit_and_figures(tmp_path, repo_root) -> None:
    study = PharmacologyStudy.from_files(EXAMPLE_DIR, mode="strict")
    destination = tmp_path / "analysis"
    manifest = study.export_analysis_bundle(
        destination,
        run_id="iv-analysis",
        nca_config=repo_root / "configs" / "nca_single_dose_iv.yaml",
        mode="strict",
    )
    assert manifest.analysis_executed is True
    assert manifest.quality_status == "research_only"
    assert "nca_parameters" in manifest.data
    assert "model_parameters" in manifest.data
    assert "predictions" in manifest.data
    assert "diagnostics" in manifest.data
    assert "concentration_time_figure" in manifest.data
    omitted = {item.role for item in manifest.omitted}
    assert "nca_parameters" not in omitted
    assert "concentration_time_figure" not in omitted
    assert "environment_lock" in omitted
    svg = (destination / "figures" / "concentration_time.svg").read_text(encoding="utf-8")
    assert "data-role='concentration-time'" in svg or 'data-role="concentration-time"' in svg
    assert "data-status='BLQ' data-time='24' data-conc=''" in svg
    assert "data-status='quantified' data-time='24'" not in svg
    assert "SUBJ-01" in svg and "SUBJ-02" in svg
    assert (destination / "figures" / "goodness_of_fit.svg").is_file()
    assert (destination / "figures" / "residual_diagnostics.svg").is_file()
    report = (destination / "reports" / "analysis_report.html").read_text(encoding="utf-8")
    assert "Research use only" in report
    assert "not a dosing recommendation" in report
    validation = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "were executed" in validation
    assert "were not run" not in validation
    nca = pq.read_table(destination / "analysis" / "nca_parameters.parquet")
    assert nca.num_rows == 2
    parameters = json.loads((destination / "provenance" / "analysis_parameters.json").read_text(encoding="utf-8"))
    assert parameters["analysis_executed"] is True
    blq = [item for item in study.observations if item.sample_id.endswith("T1440")]
    assert all(item.numeric_value is None for item in blq)


def test_failed_fit_omits_estimates_and_still_marks_the_analysis_as_run(tmp_path) -> None:
    study = build_study()
    records = AnalysisRecords(
        model_fits=[
            ModelFit(
                fit_id="FIT-SUBJ-01",
                model_id="one_compartment_iv_bolus",
                dataset_id=study.metadata.study_id,
                software_version=__version__,
                estimation_method="nonlinear least squares",
                context_of_use="Research-use fit that did not converge on this fixture.",
                diagnostic_status=DiagnosticStatus.FAIL,
                converged=False,
                limitations=["The optimizer did not converge. Estimates were omitted."],
            )
        ],
        diagnostics=[
            DiagnosticRecord(
                fit_id="FIT-SUBJ-01",
                status=DiagnosticStatus.FAIL,
                limitations=["The optimizer did not converge. Estimates were omitted."],
                warnings=["non-convergence"],
            )
        ],
        nca_parameters=[
            NcaParameterSet(
                subject_id="SUBJ-01",
                analyte_id="PW-0001-PARENT",
                route=Route.IV_BOLUS,
                dose_event_id="DOSE-SUBJ-01",
                auc_method="linear_up_log_down",
                terminal_phase_qualified=False,
                auc_last=1.0,
                auc_last_unit="mg*h/L",
                cmax=1.0,
                cmax_unit="mg/L",
                tmax=1.0,
                tmax_unit="h",
                warnings=["Terminal parameters were omitted."],
            )
        ],
    )
    payload = AnalysisPayload(
        records=records,
        concentration_time_svg="<svg xmlns='http://www.w3.org/2000/svg' data-role='concentration-time'></svg>",
        analysis_parameters={"research_only": True},
    )
    destination = tmp_path / "failed-fit"
    manifest = export_bundle(study, destination, run_id="failed-fit", analysis=payload)
    assert manifest.analysis_executed is True
    assert "nca_parameters" in manifest.data
    assert "diagnostics" in manifest.data
    omitted = {item.role: item.reason for item in manifest.omitted}
    assert "model_parameters" in omitted
    assert "omitted" in omitted["model_parameters"]
    assert "predictions" in omitted
    assert "goodness_of_fit_figure" in omitted
    assert not (destination / "analysis" / "model_parameters.parquet").exists()


def test_cli_analyze_writes_a_bundle(tmp_path, repo_root) -> None:
    bundle = tmp_path / "cli-analysis"
    status = main(
        [
            "analyze",
            str(EXAMPLE_DIR),
            "--config",
            str(repo_root / "configs" / "nca_single_dose_iv.yaml"),
            "--bundle",
            str(bundle),
            "--run-id",
            "cli-analysis",
            "--strict",
        ]
    )
    assert status == 0
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["analysis_executed"] is True
    assert manifest["quality_status"] == "research_only"
