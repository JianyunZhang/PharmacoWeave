"""Run the locked Phase 1 analysis and write its result bundle.

The public steps remain separate: ``run_nca`` calculates NCA, and
``fit_one_compartment_iv_bolus`` estimates CL and V. This module is the path
that draws figures and sets ``analysis_executed``.
"""

from pathlib import Path

from pharmacoweave.core.analysis_export import AnalysisPayload
from pharmacoweave.core.modeling import AnalysisRecords
from pharmacoweave.core.provenance import ResultBundleManifest
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
from pharmacoweave.diagnostics.plots import (
    PlotPoint,
    PlotSeries,
    concentration_time_svg,
    goodness_of_fit_svg,
    residual_svg,
)
from pharmacoweave.errors import ValidationFailed
from pharmacoweave.pk_models.fit import IvBolusFit, fit_one_compartment_iv_bolus
from pharmacoweave.pk_nca.config import NcaConfig
from pharmacoweave.pk_nca.runner import load_nca_config, run_nca
from pharmacoweave.pk_profiles import load_iv_bolus_profiles


def export_analysis_bundle(
    study: PharmacologyStudy,
    path: str | Path,
    *,
    run_id: str,
    nca_config: NcaConfig | str | Path,
    mode: str = "research",
    error_model: str = "proportional",
) -> ResultBundleManifest:
    """Validate, run NCA and the IV bolus fit, and write a new analysis bundle."""

    report = validate_study(study, mode=mode)
    if not report.ok:
        raise ValidationFailed(report)
    loaded = nca_config if isinstance(nca_config, NcaConfig) else load_nca_config(nca_config)
    nca_records = run_nca(study, loaded)
    fitted = fit_one_compartment_iv_bolus(study, error_model=error_model)
    payload = _payload(study, loaded, nca_records, fitted, error_model=error_model)
    from pharmacoweave.core.bundle import export_bundle

    return export_bundle(study, path, run_id=run_id, mode=mode, analysis=payload)


def _payload(
    study: PharmacologyStudy,
    config: NcaConfig,
    nca_records: AnalysisRecords,
    fitted: IvBolusFit,
    *,
    error_model: str,
) -> AnalysisPayload:
    profiles = load_iv_bolus_profiles(study)
    curves = {
        (curve.subject_id, curve.analyte_id, curve.dose_event_id): curve.points
        for curve in fitted.curves
    }
    series = []
    for profile in profiles:
        points = tuple(
            PlotPoint(time_h=point.time_h, status=point.status.value, concentration=point.value)
            for point in profile.points
        )
        series.append(
            PlotSeries(
                name=f"{profile.subject_id} {profile.analyte_id}",
                unit=profile.concentration_unit or "",
                lloq=profile.lloq,
                points=points,
                prediction=curves.get((profile.subject_id, profile.analyte_id, profile.dose_event_id), ()),
            )
        )
    title = f"{study.metadata.study_id} concentration-time"
    concentration = concentration_time_svg(series, title=title)
    goodness = None
    residuals = None
    if fitted.predictions:
        units = sorted({row.unit for row in fitted.predictions})
        unit = units[0] if len(units) == 1 else "reported unit"
        goodness = goodness_of_fit_svg(
            [(row.observed, row.predicted) for row in fitted.predictions],
            unit=unit,
            title=f"{study.metadata.study_id} goodness of fit",
        )
        residuals = residual_svg(
            [(row.time_h, row.weighted_residual) for row in fitted.predictions],
            title=f"{study.metadata.study_id} weighted residuals",
        )
    merged = AnalysisRecords(
        model_specs=list(fitted.records.model_specs),
        model_fits=list(fitted.records.model_fits),
        parameter_estimates=list(fitted.records.parameter_estimates),
        nca_parameters=list(nca_records.nca_parameters),
        diagnostics=list(fitted.records.diagnostics),
    )
    return AnalysisPayload(
        records=merged,
        predictions=fitted.predictions,
        concentration_time_svg=concentration,
        goodness_of_fit_svg=goodness,
        residual_diagnostics_svg=residuals,
        analysis_parameters={
            "nca": config.model_dump(mode="json"),
            "fit": {
                "model_id": "one_compartment_iv_bolus",
                "error_model": error_model,
                "estimation_method": "nonlinear least squares",
            },
            "research_only": True,
        },
    )
