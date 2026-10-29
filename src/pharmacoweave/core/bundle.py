"""Write an immutable, research-use result bundle."""

import hashlib
import json
import platform
import shutil
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path

import jsonschema
import pyarrow as pa
import pyarrow.parquet as pq

from pharmacoweave.core.analysis_export import AnalysisPayload, PredictionRow
from pharmacoweave.core.entities import ObservationType
from pharmacoweave.core.evidence import EvidenceExchange, ReviewStatus
from pharmacoweave.core.modeling import AnalysisRecords, NcaParameterSet, ParameterEstimate
from pharmacoweave.core.provenance import (
    BundleFile,
    OmittedArtifact,
    Producer,
    ResultBundleManifest,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import ValidationReport, validate_study
from pharmacoweave.errors import IngestError, ValidationFailed
from pharmacoweave.ingestion.tabular import (
    analyte_row,
    assay_row,
    compound_row,
    dose_row,
    identity_row,
    observation_row,
    sample_row,
    subject_row,
    write_named_table,
)
from pharmacoweave.reporting.html import render_validation_report
from pharmacoweave.schema_export import schema_directory
from pharmacoweave.version import SCHEMA_VERSION, __version__

_PATHS = {
    "study_metadata": "study/study_metadata.json",
    "subjects": "study/subjects.parquet",
    "compounds": "study/compounds.parquet",
    "analytes": "study/analytes.parquet",
    "assays": "study/assays.parquet",
    "identity_map": "study/identity_map.parquet",
    "dosing_events": "experimental/dosing_events.parquet",
    "biosamples": "experimental/biosamples.parquet",
    "concentration_observations": "experimental/concentration_observations.parquet",
    "pd_observations": "experimental/pd_observations.parquet",
    "toxicology_observations": "experimental/toxicology_observations.parquet",
    "source_files": "provenance/source_files.json",
    "analysis_parameters": "provenance/analysis_parameters.json",
    "software_versions": "provenance/software_versions.json",
    "runtime": "provenance/runtime.json",
    "pharmacology_evidence": "evidence/pharmacology_evidence.json",
    "cross_project_links": "evidence/cross_project_links.json",
    "ingest_validation": "qc/ingest_validation.json",
    "quality_metrics": "qc/quality_metrics.json",
    "validation_report": "reports/validation_report.html",
}

_ENVIRONMENT_OMITTED = OmittedArtifact(
    role="environment_lock",
    reason=(
        "provenance/runtime.json records the writer runtime. "
        "A sealed environment.lock is not produced."
    ),
)

_VALIDATION_OMITTED = (
    OmittedArtifact(
        role="nca_parameters",
        reason="analysis/nca_parameters.parquet is absent because NCA was not run.",
    ),
    OmittedArtifact(
        role="model_parameters",
        reason="analysis/model_parameters.parquet is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="predictions",
        reason="analysis/predictions.parquet is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="residuals",
        reason="analysis/residuals.parquet is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="diagnostics",
        reason="analysis/diagnostics.json is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="concentration_time_figure",
        reason="figures/concentration_time.svg is absent because no analysis figure was drawn.",
    ),
    OmittedArtifact(
        role="goodness_of_fit_figure",
        reason="figures/goodness_of_fit.svg is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="residual_diagnostics_figure",
        reason="figures/residual_diagnostics.svg is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="analysis_report",
        reason=(
            "reports/analysis_report.html is absent. reports/validation_report.html is the data-validation inventory."
        ),
    ),
)

_ANALYSIS_NOTE = "No NCA or model fit was run."
_LINK_POLICY = (
    "Identical subject_id text across projects is not a match. "
    "Only review_status=accepted is eligible for a later join. This bundle does not join records."
)


def export_bundle(
    study: PharmacologyStudy,
    path: str | Path,
    *,
    run_id: str,
    mode: str = "research",
    analysis: AnalysisPayload | None = None,
) -> ResultBundleManifest:
    """Validate ``study`` and write a new bundle directory.

    Without ``analysis`` the bundle is a validation inventory and
    ``analysis_executed`` stays false. With ``analysis``, calculated artifacts
    are written and ``analysis_executed`` is true.
    """

    report = validate_study(study, mode=mode)
    if not report.ok:
        raise ValidationFailed(report)
    destination = Path(path)
    if destination.exists():
        raise IngestError(f"{destination} already exists. Each run writes a new directory.")
    destination.mkdir(parents=True)
    try:
        return _write_bundle(study, destination, run_id=run_id, report=report, analysis=analysis)
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _write_bundle(
    study: PharmacologyStudy,
    destination: Path,
    *,
    run_id: str,
    report: ValidationReport,
    analysis: AnalysisPayload | None,
) -> ResultBundleManifest:
    created_at = datetime.now(UTC)
    counts = _counts(study, report)
    paths = dict(_PATHS)
    _write_study_tables(study, destination)
    if analysis is None:
        omitted = [_ENVIRONMENT_OMITTED, *_VALIDATION_OMITTED]
        parameters_doc: dict[str, object] = {
            "analysis_executed": False,
            "parameters": {},
            "note": _ANALYSIS_NOTE,
        }
        analysis_executed = False
    else:
        omitted = [_ENVIRONMENT_OMITTED]
        _stage_analysis(destination, analysis, paths, omitted)
        parameters_doc = {
            "analysis_executed": True,
            "parameters": analysis.analysis_parameters,
            "note": (
                "NCA and the one-compartment IV bolus fit were executed when this bundle was written. "
                "Omitted artifacts are listed in the manifest. Research use only."
            ),
        }
        analysis_executed = True
        report_html = _analysis_report(
            study,
            analysis.records,
            run_id=run_id,
            generated_at=created_at,
            omitted=omitted,
        )
        relative = "reports/analysis_report.html"
        _write_text(destination / relative, report_html)
        paths["analysis_report"] = relative

    _write_json(
        destination / paths["source_files"],
        {"files": [item.model_dump(mode="json") for item in study.source_files]},
    )
    _write_json(destination / paths["analysis_parameters"], parameters_doc)
    _write_json(destination / paths["software_versions"], _software_versions())
    _write_json(
        destination / paths["runtime"],
        {
            "python": sys.version,
            "platform": platform.platform(),
            "note": "Runtime that wrote this bundle. Not a sealed analysis environment.",
        },
    )
    exchange = EvidenceExchange(records=list(study.external_evidence))
    _write_json(destination / paths["pharmacology_evidence"], json.loads(exchange.model_dump_json()))
    _write_json(
        destination / paths["cross_project_links"],
        {
            "policy": _LINK_POLICY,
            "links": [item.model_dump(mode="json") for item in study.identity_links],
        },
    )
    validation_payload = report.model_dump(mode="json")
    validation_payload["ok"] = report.ok
    _write_json(destination / paths["ingest_validation"], validation_payload)
    _write_json(destination / paths["quality_metrics"], counts)
    html = render_validation_report(
        study,
        report,
        run_id=run_id,
        generated_at=created_at,
        counts=counts,
        omitted=omitted,
        version=__version__,
        schema_version=SCHEMA_VERSION,
        analysis_executed=analysis_executed,
    )
    _write_text(destination / paths["validation_report"], html)

    files = [
        BundleFile(
            path=relative,
            role=role,
            sha256=_sha256(destination / relative),
            media_type=_media_type(relative),
        )
        for role, relative in paths.items()
    ]
    manifest = ResultBundleManifest(
        run_id=run_id,
        study_id=study.metadata.study_id,
        created_at=created_at,
        producer=Producer(version=__version__),
        files=files,
        data=paths,
        omitted=omitted,
        analysis_executed=analysis_executed,
        warnings=[issue.message for issue in report.issues if issue.severity == "warning"],
    )
    payload = json.loads(manifest.model_dump_json())
    _check_manifest_schema(payload)
    manifest_path = destination / "manifest.json"
    manifest_bytes = (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    manifest_path.write_bytes(manifest_bytes)
    digest = hashlib.sha256(manifest_bytes).hexdigest()
    (destination / "manifest.sha256").write_text(f"{digest}  manifest.json\n", encoding="utf-8")
    return manifest


def _write_study_tables(study: PharmacologyStudy, destination: Path) -> None:
    _write_json(destination / _PATHS["study_metadata"], json.loads(study.metadata.model_dump_json()))
    write_named_table(destination / _PATHS["subjects"], "subjects", [subject_row(item) for item in study.subjects])
    write_named_table(destination / _PATHS["compounds"], "compounds", [compound_row(item) for item in study.compounds])
    write_named_table(destination / _PATHS["analytes"], "analytes", [analyte_row(item) for item in study.analytes])
    write_named_table(destination / _PATHS["assays"], "assays", [assay_row(item) for item in study.assays])
    write_named_table(
        destination / _PATHS["identity_map"],
        "identity_map",
        [identity_row(item) for item in study.identity_links],
    )
    write_named_table(
        destination / _PATHS["dosing_events"],
        "dosing_events",
        [dose_row(item) for item in study.dose_events],
    )
    write_named_table(
        destination / _PATHS["biosamples"],
        "biosamples",
        [sample_row(item) for item in study.biosamples],
    )
    by_type = {
        ObservationType.CONCENTRATION: "concentration_observations",
        ObservationType.EFFECT: "pd_observations",
        ObservationType.TOXICITY: "toxicology_observations",
    }
    partitioned: dict[ObservationType, list[dict[str, object]]] = {kind: [] for kind in by_type}
    for observation in study.observations:
        partitioned[observation.observation_type].append(observation_row(observation))
    if sum(len(rows) for rows in partitioned.values()) != len(study.observations):
        raise IngestError("Every observation must be written to exactly one experimental table.")
    for kind, role in by_type.items():
        write_named_table(destination / _PATHS[role], "observations", partitioned[kind])


def _stage_analysis(
    destination: Path,
    analysis: AnalysisPayload,
    paths: dict[str, str],
    omitted: list[OmittedArtifact],
) -> None:
    records = analysis.records
    _check_analysis_schema(records)
    relative = "analysis/analysis_records.json"
    _write_json(destination / relative, json.loads(records.model_dump_json()))
    paths["analysis_records"] = relative
    if records.nca_parameters:
        relative = "analysis/nca_parameters.parquet"
        _write_rows(destination / relative, [_nca_row(item) for item in records.nca_parameters], _NCA_TYPES)
        paths["nca_parameters"] = relative
    else:
        omitted.append(
            OmittedArtifact(
                role="nca_parameters",
                reason="analysis/nca_parameters.parquet is absent because NCA returned no profiles.",
            )
        )
    if records.parameter_estimates:
        relative = "analysis/model_parameters.parquet"
        _write_rows(
            destination / relative,
            [_estimate_row(item) for item in records.parameter_estimates],
            _ESTIMATE_TYPES,
        )
        paths["model_parameters"] = relative
    else:
        omitted.append(
            OmittedArtifact(
                role="model_parameters",
                reason=(
                    "analysis/model_parameters.parquet is absent because the fit failed, "
                    "was not identifiable, or hit a bound. Estimates were omitted."
                ),
            )
        )
    if analysis.predictions:
        relative = "analysis/predictions.parquet"
        _write_rows(destination / relative, [_prediction_row(item) for item in analysis.predictions], _PREDICTION_TYPES)
        paths["predictions"] = relative
        relative = "analysis/residuals.parquet"
        _write_rows(destination / relative, [_residual_row(item) for item in analysis.predictions], _RESIDUAL_TYPES)
        paths["residuals"] = relative
    else:
        omitted.append(
            OmittedArtifact(
                role="predictions",
                reason="analysis/predictions.parquet is absent because no accepted fit produced predictions.",
            )
        )
        omitted.append(
            OmittedArtifact(
                role="residuals",
                reason="analysis/residuals.parquet is absent because no accepted fit produced residuals.",
            )
        )
    if records.diagnostics or records.model_fits:
        relative = "analysis/diagnostics.json"
        _write_json(
            destination / relative,
            {
                "research_only": True,
                "note": "A passing diagnostic states limitations. A failed fit does not carry parameter estimates.",
                "model_fits": [item.model_dump(mode="json") for item in records.model_fits],
                "diagnostics": [item.model_dump(mode="json") for item in records.diagnostics],
            },
        )
        paths["diagnostics"] = relative
    else:
        omitted.append(
            OmittedArtifact(
                role="diagnostics",
                reason="analysis/diagnostics.json is absent because this analysis did not record a fit diagnostic.",
            )
        )
    _stage_svg(
        destination,
        paths,
        omitted,
        role="concentration_time_figure",
        relative="figures/concentration_time.svg",
        svg=analysis.concentration_time_svg,
        absent="figures/concentration_time.svg is absent because no concentration-time figure was produced.",
    )
    _stage_svg(
        destination,
        paths,
        omitted,
        role="goodness_of_fit_figure",
        relative="figures/goodness_of_fit.svg",
        svg=analysis.goodness_of_fit_svg or "",
        absent="figures/goodness_of_fit.svg is absent because no accepted fit was produced.",
    )
    _stage_svg(
        destination,
        paths,
        omitted,
        role="residual_diagnostics_figure",
        relative="figures/residual_diagnostics.svg",
        svg=analysis.residual_diagnostics_svg or "",
        absent="figures/residual_diagnostics.svg is absent because no accepted fit was produced.",
    )


def _stage_svg(
    destination: Path,
    paths: dict[str, str],
    omitted: list[OmittedArtifact],
    *,
    role: str,
    relative: str,
    svg: str,
    absent: str,
) -> None:
    if svg.strip():
        _write_text(destination / relative, svg)
        paths[role] = relative
        return
    omitted.append(OmittedArtifact(role=role, reason=absent))


def _analysis_report(
    study: PharmacologyStudy,
    records: AnalysisRecords,
    *,
    run_id: str,
    generated_at: datetime,
    omitted: list[OmittedArtifact],
) -> str:
    from pharmacoweave.reporting.analysis_html import render_analysis_report

    return render_analysis_report(
        records,
        study_id=study.metadata.study_id,
        title=study.metadata.title,
        run_id=run_id,
        generated_at=generated_at,
        omitted=omitted,
        version=__version__,
        schema_version=SCHEMA_VERSION,
    )


def _counts(study: PharmacologyStudy, report: ValidationReport) -> dict[str, int]:
    observations = study.observations
    issues = report.issues
    return {
        "n_subjects": len(study.subjects),
        "n_dose_events": len(study.dose_events),
        "n_biosamples": len(study.biosamples),
        "n_observations": len(observations),
        "n_quantified": sum(item.status.value == "quantified" for item in observations),
        "n_blq": sum(item.status.value == "BLQ" for item in observations),
        "n_above_uloq": sum(item.status.value == "above_ULOQ" for item in observations),
        "n_missing": sum(item.status.value == "missing" for item in observations),
        "n_invalid": sum(item.status.value == "invalid" for item in observations),
        "n_identity_links_accepted": sum(item.review_status is ReviewStatus.ACCEPTED for item in study.identity_links),
        "n_external_evidence": len(study.external_evidence),
        "n_errors": sum(item.severity == "error" for item in issues),
        "n_warnings": sum(item.severity == "warning" for item in issues),
    }


def _software_versions() -> dict[str, str]:
    return {
        "pharmacoweave": __version__,
        "python": platform.python_version(),
        "pydantic": metadata.version("pydantic"),
        "pyarrow": metadata.version("pyarrow"),
        "pandas": metadata.version("pandas"),
        "jinja2": metadata.version("jinja2"),
        "numpy": metadata.version("numpy"),
        "scipy": metadata.version("scipy"),
    }


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _check_manifest_schema(payload: dict[str, object]) -> None:
    schema_path = schema_directory() / "result_bundle.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    jsonschema.validate(instance=payload, schema=schema, format_checker=jsonschema.FormatChecker())


def _check_analysis_schema(records: AnalysisRecords) -> None:
    schema_path = schema_directory() / "analysis_records.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    jsonschema.validate(instance=json.loads(records.model_dump_json()), schema=schema)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _media_type(path: str) -> str:
    if path.endswith(".json"):
        return "application/json"
    if path.endswith(".parquet"):
        return "application/vnd.apache.parquet"
    if path.endswith(".html"):
        return "text/html"
    if path.endswith(".svg"):
        return "image/svg+xml"
    return "application/octet-stream"


def _nca_row(record: NcaParameterSet) -> dict[str, object]:
    payload = record.model_dump(mode="json")
    payload["warnings"] = "|".join(part.replace("|", "/") for part in record.warnings)
    return payload


def _estimate_row(record: ParameterEstimate) -> dict[str, object]:
    return record.model_dump(mode="json")


def _prediction_row(row: PredictionRow) -> dict[str, object]:
    return {
        "fit_id": row.fit_id,
        "observation_id": row.observation_id,
        "subject_id": row.subject_id,
        "time_h": row.time_h,
        "observed": row.observed,
        "predicted": row.predicted,
        "residual": row.residual,
        "weighted_residual": row.weighted_residual,
        "unit": row.unit,
    }


def _residual_row(row: PredictionRow) -> dict[str, object]:
    return {
        "fit_id": row.fit_id,
        "observation_id": row.observation_id,
        "subject_id": row.subject_id,
        "time_h": row.time_h,
        "residual": row.residual,
        "weighted_residual": row.weighted_residual,
        "unit": row.unit,
    }


def _write_rows(path: Path, rows: list[dict[str, object]], types: tuple[tuple[str, pa.DataType], ...]) -> None:
    arrays = [pa.array([row.get(name) for row in rows], type=datatype) for name, datatype in types]
    table = pa.Table.from_arrays(arrays, names=[name for name, _datatype in types])
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path)


_STRING = pa.string()
_FLOAT = pa.float64()
_BOOL = pa.bool_()
_NCA_TYPES = (
    ("subject_id", _STRING),
    ("analyte_id", _STRING),
    ("route", _STRING),
    ("dose_event_id", _STRING),
    ("auc_method", _STRING),
    ("terminal_phase_qualified", _BOOL),
    ("auc_last", _FLOAT),
    ("auc_last_unit", _STRING),
    ("auc_inf", _FLOAT),
    ("auc_inf_unit", _STRING),
    ("auc_inf_extrap_fraction", _FLOAT),
    ("cmax", _FLOAT),
    ("cmax_unit", _STRING),
    ("tmax", _FLOAT),
    ("tmax_unit", _STRING),
    ("lambda_z", _FLOAT),
    ("lambda_z_unit", _STRING),
    ("half_life", _FLOAT),
    ("half_life_unit", _STRING),
    ("terminal_phase_n_points", pa.int64()),
    ("terminal_phase_adjusted_r2", _FLOAT),
    ("clearance_semantic", _STRING),
    ("clearance", _FLOAT),
    ("clearance_unit", _STRING),
    ("volume_semantic", _STRING),
    ("volume", _FLOAT),
    ("volume_unit", _STRING),
    ("excluded_reason", _STRING),
    ("warnings", _STRING),
)
_ESTIMATE_TYPES = (
    ("fit_id", _STRING),
    ("parameter", _STRING),
    ("semantic", _STRING),
    ("value", _FLOAT),
    ("unit", _STRING),
    ("estimation_method", _STRING),
    ("route", _STRING),
    ("standard_error", _FLOAT),
    ("ci_lower", _FLOAT),
    ("ci_upper", _FLOAT),
    ("ci_level", _FLOAT),
    ("bound_hit", _BOOL),
)
_PREDICTION_TYPES = (
    ("fit_id", _STRING),
    ("observation_id", _STRING),
    ("subject_id", _STRING),
    ("time_h", _FLOAT),
    ("observed", _FLOAT),
    ("predicted", _FLOAT),
    ("residual", _FLOAT),
    ("weighted_residual", _FLOAT),
    ("unit", _STRING),
)
_RESIDUAL_TYPES = (
    ("fit_id", _STRING),
    ("observation_id", _STRING),
    ("subject_id", _STRING),
    ("time_h", _FLOAT),
    ("residual", _FLOAT),
    ("weighted_residual", _FLOAT),
    ("unit", _STRING),
)
