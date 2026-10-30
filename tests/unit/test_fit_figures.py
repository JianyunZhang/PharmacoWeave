"""Diagnostic SVG figures omit BLQ and are stored when analyses are exported."""

import json
import xml.etree.ElementTree as ET

import pyarrow.parquet as pq
from examples.direct_effect_and_evidence import E0, build_effects, build_exposures
from examples.iv_bolus_case import build_study

from pharmacoweave.core.profiles import iv_bolus_profiles
from pharmacoweave.pk_models.fit_iv import fit_iv_bolus
from pharmacoweave.pk_models.wnls import fit_iv_bolus_wnls
from pharmacoweave.pkpd import fit_linked_direct_emax
from pharmacoweave.reporting.figures import (
    EffectCurveOverlay,
    EffectPlotPoint,
    PredictedPoint,
    SemilogOverlay,
    concentration_effect_svg,
    concentration_time_semilog_svg,
    goodness_of_fit_svg,
    residual_diagnostics_svg,
    weighted_residual_svg,
)


def test_non_positive_observations_are_not_plotted() -> None:
    points = [
        PredictedPoint("SUBJ-01", 1.0, 1.0, 0.9),
        PredictedPoint("SUBJ-01", 24.0, 0.0, 0.01),
    ]
    for figure in (
        goodness_of_fit_svg(points),
        residual_diagnostics_svg(points),
        weighted_residual_svg(points),
    ):
        ET.fromstring(figure)
        assert figure.count("<circle") == 1
        assert "not plotted as zero" in figure


def test_semilog_and_concentration_effect_svgs_are_valid_xml_without_blq() -> None:
    study = build_study()
    profiles = iv_bolus_profiles(study)
    fits = fit_iv_bolus(study)
    overlays = [
        SemilogOverlay(item.subject_id, item.times_h, item.predicted_mg_l)
        for item in fits
        if item.predicted_mg_l
    ]
    semilog = concentration_time_semilog_svg(
        profiles,
        lloq_mg_l=0.01,
        overlays=overlays,
        terminal_point_counts={"SUBJ-01": 4, "SUBJ-02": 4},
    )
    ET.fromstring(semilog)
    assert "LLOQ" in semilog
    assert "Semi-log concentration vs time" in semilog
    assert semilog.count('class="observed-point"') == 16
    assert semilog.count('class="terminal-point"') == 8

    _linked, emax_fit = fit_linked_direct_emax(build_effects(), build_exposures(), e0=E0, e0_unit="percent")
    eff_points = [
        EffectPlotPoint("direct-emax", pair.concentration, pair.effect) for pair in emax_fit.pairs
    ]
    eff_Points_with_zero = [*eff_points, EffectPlotPoint("direct-emax", 0.0, 0.0)]
    eff_svg = concentration_effect_svg(
        eff_Points_with_zero,
        overlays=[EffectCurveOverlay("direct-emax", e0=10.0, emax=80.0, ec50=8.0, hill=1.0)],
        concentration_unit="mg/L",
        effect_unit="percent",
    )
    ET.fromstring(eff_svg)
    assert eff_svg.count("<circle") == len(eff_points)
    assert "EC50" in eff_svg
    assert 'class="emax-curve"' in eff_svg


def test_fit_bundle_writes_figures_and_uncertainty_without_blq(tmp_path) -> None:
    study = build_study()
    fits = fit_iv_bolus_wnls(study, weight="1/y^2")
    destination = tmp_path / "fit-bundle"
    manifest = study.export_bundle(destination, run_id="fit-run", fits=fits)
    assert manifest.analysis_executed is True
    assert manifest.data["goodness_of_fit_figure"] == "figures/goodness_of_fit.svg"
    assert manifest.data["residual_diagnostics_figure"] == "figures/residual_diagnostics.svg"
    assert manifest.data["weighted_residuals_figure"] == "figures/weighted_residuals.svg"
    assert manifest.data["concentration_time_semilog_figure"] == "figures/concentration_time_semilog.svg"
    assert not any(item.role == "goodness_of_fit_figure" for item in manifest.omitted)

    gof = (destination / "figures" / "goodness_of_fit.svg").read_text(encoding="utf-8")
    residual = (destination / "figures" / "residual_diagnostics.svg").read_text(encoding="utf-8")
    wres = (destination / "figures" / "weighted_residuals.svg").read_text(encoding="utf-8")
    assert gof.count("<circle") == 16
    assert residual.count("<circle") == 16
    assert wres.count("<circle") == 16
    assert "Observed vs predicted" in gof
    assert "Residual vs time" in residual
    assert "Weighted residual vs time" in wres
    assert "not plotted as zero" in gof

    params = pq.read_table(destination / "analysis" / "model_parameters.parquet")
    assert "standard_error" in params.column_names
    assert "ci_lower" in params.column_names
    assert "ci_upper" in params.column_names
    assert all(val == 0.95 for val in params.column("ci_level").to_pylist())

    predictions = pq.read_table(destination / "analysis" / "predictions.parquet")
    assert predictions.num_rows == 16
    assert min(predictions.column("observed_mg_l").to_pylist()) > 0
    assert all("T1440" not in observation_id for observation_id in predictions.column("observation_id").to_pylist())
    diagnostics = json.loads((destination / "analysis" / "diagnostics.json").read_text(encoding="utf-8"))
    assert diagnostics["research_use_only"] is True
    html = (destination / "reports" / "validation_report.html").read_text(encoding="utf-8")
    assert "A one-compartment IV bolus fit was run." in html
    assert "Observed vs predicted" in html
    assert "Weighted residual vs time" in html


def test_nca_only_bundle_does_not_invent_a_fit_figure(tmp_path, repo_root) -> None:
    study = build_study()
    results = study.run_nca(repo_root / "configs" / "nca_single_dose_iv.yaml")
    destination = tmp_path / "nca-only"
    manifest = study.export_bundle(destination, run_id="nca-only", nca_parameters=results)
    assert any(item.role == "goodness_of_fit_figure" for item in manifest.omitted)
    assert not (destination / "figures" / "goodness_of_fit.svg").exists()
    assert (destination / "figures" / "concentration_time.svg").is_file()
    assert (destination / "figures" / "concentration_time_semilog.svg").is_file()
