"""Goodness-of-fit SVG omits BLQ and is stored when a fit is exported."""

import json

import pyarrow.parquet as pq
from examples.iv_bolus_case import build_study

from pharmacoweave.pk_models.fit_iv import fit_iv_bolus
from pharmacoweave.reporting.figures import PredictedPoint, goodness_of_fit_svg, residual_diagnostics_svg


def test_non_positive_observations_are_not_plotted() -> None:
    points = [
        PredictedPoint("SUBJ-01", 1.0, 1.0, 0.9),
        PredictedPoint("SUBJ-01", 24.0, 0.0, 0.01),
    ]
    for figure in (goodness_of_fit_svg(points), residual_diagnostics_svg(points)):
        assert figure.count("<circle") == 1
        assert "not plotted as zero" in figure
        assert "Observed" in figure or "Residual" in figure


def test_fit_bundle_writes_figures_without_blq(tmp_path) -> None:
    study = build_study()
    fits = fit_iv_bolus(study)
    destination = tmp_path / "fit-bundle"
    manifest = study.export_bundle(destination, run_id="fit-run", fits=fits)
    assert manifest.analysis_executed is True
    assert manifest.data["goodness_of_fit_figure"] == "figures/goodness_of_fit.svg"
    assert manifest.data["residual_diagnostics_figure"] == "figures/residual_diagnostics.svg"
    assert not any(item.role == "goodness_of_fit_figure" for item in manifest.omitted)
    gof = (destination / "figures" / "goodness_of_fit.svg").read_text(encoding="utf-8")
    residual = (destination / "figures" / "residual_diagnostics.svg").read_text(encoding="utf-8")
    assert gof.count("<circle") == 16
    assert residual.count("<circle") == 16
    assert "Observed vs predicted" in gof
    assert "Residual vs time" in residual
    assert "not plotted as zero" in gof
    predictions = pq.read_table(destination / "analysis" / "predictions.parquet")
    assert predictions.num_rows == 16
    assert min(predictions.column("observed_mg_l").to_pylist()) > 0
    assert all("T1440" not in observation_id for observation_id in predictions.column("observation_id").to_pylist())
    diagnostics = json.loads((destination / "analysis" / "diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["research_use_only"] is True
    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "A one-compartment IV bolus fit was run." in html
    assert "Observed vs predicted" in html


def test_nca_only_bundle_does_not_invent_a_fit_figure(tmp_path, repo_root) -> None:
    study = build_study()
    results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    destination = tmp_path / "nca-only"
    manifest = study.export_bundle(destination, run_id="nca-only", nca_parameters=results)
    assert any(item.role == "goodness_of_fit_figure" for item in manifest.omitted)
    assert not (destination / "figures" / "goodness_of_fit.svg").exists()
