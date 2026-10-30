"""Write an immutable, research-use result bundle."""

import hashlib
import json
import platform
import shutil
import sys
from collections.abc import Sequence
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

import jsonschema

from pharmacoweave.core.entities import ObservationType
from pharmacoweave.core.evidence import EvidenceExchange, ReviewStatus
from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.core.provenance import (
    BundleFile,
    OmittedArtifact,
    Producer,
    ResultBundleManifest,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.validation import validate_study
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
from pharmacoweave.reporting.figures import (
    PredictedPoint,
    concentration_time_svg,
    goodness_of_fit_svg,
    residual_diagnostics_svg,
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

_OMITTED = (
    OmittedArtifact(
        role="environment_lock",
        reason=(
            "provenance/runtime.json records the writer runtime. "
            "A sealed environment.lock is not produced for a validation-only bundle."
        ),
    ),
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
    nca_parameters: list[NcaParameterSet] | None = None,
    fits: Sequence[Any] | None = None,
) -> ResultBundleManifest:
    """Validate ``study`` and write a new bundle directory.

    Pass ``nca_parameters`` from ``run_nca`` to store NCA. Pass IV bolus
    ``fits`` to store parameters, residuals, and goodness-of-fit figures.
    Omitting both leaves ``analysis_executed`` false. BLQ points are not
    drawn; a non-positive concentration is left out of the figures.
    """

    report = validate_study(study, mode=mode)
    if not report.ok:
        raise ValidationFailed(report)
    destination = Path(path)
    if destination.exists():
        raise IngestError(f"{destination} already exists. Each run writes a new directory.")
    destination.mkdir(parents=True)
    try:
        return _write_bundle(
            study,
            destination,
            run_id=run_id,
            report=report,
            nca_parameters=nca_parameters,
            fits=fits,
        )
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _write_bundle(
    study,
    destination: Path,
    *,
    run_id: str,
    report,
    nca_parameters: list[NcaParameterSet] | None,
    fits: Sequence[Any] | None,
) -> ResultBundleManifest:
    created_at = datetime.now(UTC)
    counts = _counts(study, report)
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

    _write_json(
        destination / _PATHS["source_files"],
        {"files": [item.model_dump(mode="json") for item in study.source_files]},
    )
    analysis_executed = nca_parameters is not None or fits is not None
    notes: list[str] = []
    if nca_parameters is not None:
        notes.append("Single-dose IV bolus NCA was run.")
    if fits is not None:
        notes.append("One-compartment IV bolus fit was run.")
    _write_json(
        destination / _PATHS["analysis_parameters"],
        {
            "analysis_executed": analysis_executed,
            "parameters": {},
            "note": " ".join(notes) if notes else _ANALYSIS_NOTE,
            "terminal_phase_rule": None if nca_parameters is None else "pharmacoweave.lambda_z.adjusted_r2.v1",
            "estimation_methods": [] if fits is None else sorted({item.fit.estimation_method for item in fits}),
        },
    )
    _write_json(destination / _PATHS["software_versions"], _software_versions())
    _write_json(
        destination / _PATHS["runtime"],
        {
            "python": sys.version,
            "platform": platform.platform(),
            "note": "Runtime that wrote this validation bundle. Not a sealed analysis environment.",
        },
    )
    exchange = EvidenceExchange(records=list(study.external_evidence))
    _write_json(destination / _PATHS["pharmacology_evidence"], json.loads(exchange.model_dump_json()))
    _write_json(
        destination / _PATHS["cross_project_links"],
        {
            "policy": _LINK_POLICY,
            "links": [item.model_dump(mode="json") for item in study.identity_links],
        },
    )
    validation_payload = report.model_dump(mode="json")
    validation_payload["ok"] = report.ok
    _write_json(destination / _PATHS["ingest_validation"], validation_payload)
    _write_json(destination / _PATHS["quality_metrics"], counts)
    paths = dict(_PATHS)
    omitted = list(_OMITTED)
    fit_points: list[PredictedPoint] = []
    if nca_parameters is not None:
        from pharmacoweave.core.profiles import iv_bolus_profiles

        nca_path = "analysis/nca_parameters.parquet"
        figure_path = "figures/concentration_time.svg"
        paths["nca_parameters"] = nca_path
        paths["concentration_time_figure"] = figure_path
        _write_nca_parquet(destination / nca_path, nca_parameters)
        figure = destination / figure_path
        figure.parent.mkdir(parents=True, exist_ok=True)
        figure.write_text(concentration_time_svg(iv_bolus_profiles(study)), encoding="utf-8")
        omitted = [item for item in omitted if item.role not in {"nca_parameters", "concentration_time_figure"}]
    if fits is not None:
        fit_points = _predicted_points(fits)
        fit_paths = {
            "model_parameters": "analysis/model_parameters.parquet",
            "predictions": "analysis/predictions.parquet",
            "residuals": "analysis/residuals.parquet",
            "diagnostics": "analysis/diagnostics.json",
        }
        paths.update(fit_paths)
        _write_fit_tables(destination, fits)
        omitted = [item for item in omitted if item.role not in fit_paths]
        if fit_points:
            gof_path = "figures/goodness_of_fit.svg"
            residual_path = "figures/residual_diagnostics.svg"
            paths["goodness_of_fit_figure"] = gof_path
            paths["residual_diagnostics_figure"] = residual_path
            gof_file = destination / gof_path
            residual_file = destination / residual_path
            gof_file.parent.mkdir(parents=True, exist_ok=True)
            gof_file.write_text(goodness_of_fit_svg(fit_points), encoding="utf-8")
            residual_file.write_text(residual_diagnostics_svg(fit_points), encoding="utf-8")
            omitted = [
                item
                for item in omitted
                if item.role not in {"goodness_of_fit_figure", "residual_diagnostics_figure"}
            ]
        else:
            omitted = [_figure_omitted(item) for item in omitted]
    html = render_validation_report(
        study,
        report,
        run_id=run_id,
        generated_at=created_at,
        counts=counts,
        omitted=omitted,
        version=__version__,
        schema_version=SCHEMA_VERSION,
        nca_parameters=nca_parameters,
        fit_rows=_fit_rows(fits),
        goodness_of_fit_svg=None if not fit_points else goodness_of_fit_svg(fit_points),
        residual_diagnostics_svg=None if not fit_points else residual_diagnostics_svg(fit_points),
    )
    report_path = destination / _PATHS["validation_report"]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(html, encoding="utf-8")

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


def _counts(study: PharmacologyStudy, report) -> dict[str, int]:
    observations = study.observations
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
        "n_errors": sum(item.severity == "error" for item in report.issues),
        "n_warnings": sum(item.severity == "warning" for item in report.issues),
    }


def _software_versions() -> dict[str, str]:
    return {
        "pharmacoweave": __version__,
        "python": platform.python_version(),
        "pydantic": metadata.version("pydantic"),
        "pyarrow": metadata.version("pyarrow"),
        "pandas": metadata.version("pandas"),
        "jinja2": metadata.version("jinja2"),
    }


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _check_manifest_schema(payload: dict[str, object]) -> None:
    schema_path = schema_directory() / "result_bundle.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    jsonschema.validate(instance=payload, schema=schema, format_checker=jsonschema.FormatChecker())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_nca_parquet(path: Path, records: list[NcaParameterSet]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    columns = (
        "subject_id",
        "analyte_id",
        "route",
        "dose_event_id",
        "auc_method",
        "terminal_phase_qualified",
        "auc_last",
        "auc_last_unit",
        "auc_inf",
        "auc_inf_unit",
        "auc_inf_extrap_fraction",
        "cmax",
        "cmax_unit",
        "tmax",
        "tmax_unit",
        "lambda_z",
        "lambda_z_unit",
        "half_life",
        "half_life_unit",
        "terminal_phase_n_points",
        "terminal_phase_adjusted_r2",
        "clearance_semantic",
        "clearance",
        "clearance_unit",
        "volume_semantic",
        "volume",
        "volume_unit",
        "excluded_reason",
        "warnings",
    )
    floats = {
        "auc_last",
        "auc_inf",
        "auc_inf_extrap_fraction",
        "cmax",
        "tmax",
        "lambda_z",
        "half_life",
        "terminal_phase_adjusted_r2",
        "clearance",
        "volume",
    }
    ints = {"terminal_phase_n_points"}
    rows: list[dict[str, object]] = []
    for record in records:
        payload = record.model_dump(mode="json")
        payload["warnings"] = " | ".join(record.warnings)
        rows.append({column: payload.get(column) for column in columns})
    fields = []
    for name in columns:
        if name in floats:
            datatype: pa.DataType = pa.float64()
        elif name in ints:
            datatype = pa.int64()
        elif name == "terminal_phase_qualified":
            datatype = pa.bool_()
        else:
            datatype = pa.string()
        fields.append(pa.field(name, datatype, nullable=name != "terminal_phase_qualified"))
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=pa.schema(fields)), path)


def _figure_omitted(item: OmittedArtifact) -> OmittedArtifact:
    if item.role == "goodness_of_fit_figure":
        return OmittedArtifact(
            role=item.role,
            reason="figures/goodness_of_fit.svg is absent because the fit produced no positive predictions.",
        )
    if item.role == "residual_diagnostics_figure":
        return OmittedArtifact(
            role=item.role,
            reason="figures/residual_diagnostics.svg is absent because the fit produced no positive predictions.",
        )
    return item


def _predicted_points(fits: Sequence[Any]) -> list[PredictedPoint]:
    points: list[PredictedPoint] = []
    for item in fits:
        predicted = tuple(item.predicted_mg_l)
        if not predicted:
            continue
        aligned = zip(item.times_h, item.observed_mg_l, predicted, strict=True)
        for time_h, observed, pred in aligned:
            if observed > 0 and pred > 0:
                points.append(PredictedPoint(item.subject_id, time_h, observed, pred))
    return points


def _fit_rows(fits: Sequence[Any] | None) -> list[dict[str, str]]:
    if not fits:
        return []
    rows: list[dict[str, str]] = []
    for item in fits:
        values = {estimate.parameter: estimate.value for estimate in item.estimates}
        rows.append(
            {
                "subject_id": item.subject_id,
                "method": item.fit.estimation_method,
                "status": item.diagnostic.status.value,
                "clearance": _fmt_fit(values.get("CL")),
                "volume": _fmt_fit(values.get("V")),
            }
        )
    return rows


def _fmt_fit(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def _write_fit_tables(destination: Path, fits: Sequence[Any]) -> None:
    parameter_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    residual_rows: list[dict[str, object]] = []
    for item in fits:
        for estimate in item.estimates:
            parameter_rows.append(
                {
                    "fit_id": estimate.fit_id,
                    "parameter": estimate.parameter,
                    "semantic": estimate.semantic.value,
                    "value": estimate.value,
                    "unit": estimate.unit,
                    "estimation_method": estimate.estimation_method,
                    "route": estimate.route.value,
                }
            )
        predicted = tuple(item.predicted_mg_l)
        if not predicted:
            continue
        aligned = zip(item.observation_ids, item.times_h, item.observed_mg_l, predicted, strict=True)
        for observation_id, time_h, observed, pred in aligned:
            if not (observed > 0 and pred > 0):
                continue
            residual = observed - pred
            prediction_rows.append(
                {
                    "fit_id": item.fit.fit_id,
                    "observation_id": observation_id,
                    "time_h": time_h,
                    "observed_mg_l": observed,
                    "predicted_mg_l": pred,
                    "concentration_unit": "mg/L",
                }
            )
            residual_rows.append(
                {
                    "fit_id": item.fit.fit_id,
                    "observation_id": observation_id,
                    "time_h": time_h,
                    "residual_mg_l": residual,
                    "residual_unit": "mg/L",
                }
            )
    _write_typed_parquet(
        destination / "analysis" / "model_parameters.parquet",
        (
            "fit_id",
            "parameter",
            "semantic",
            "value",
            "unit",
            "estimation_method",
            "route",
        ),
        {"value"},
        parameter_rows,
    )
    _write_typed_parquet(
        destination / "analysis" / "predictions.parquet",
        (
            "fit_id",
            "observation_id",
            "time_h",
            "observed_mg_l",
            "predicted_mg_l",
            "concentration_unit",
        ),
        {"time_h", "observed_mg_l", "predicted_mg_l"},
        prediction_rows,
    )
    _write_typed_parquet(
        destination / "analysis" / "residuals.parquet",
        ("fit_id", "observation_id", "time_h", "residual_mg_l", "residual_unit"),
        {"time_h", "residual_mg_l"},
        residual_rows,
    )
    _write_json(
        destination / "analysis" / "diagnostics.json",
        {
            "research_use_only": True,
            "note": (
                "A passing diagnostic states limitations. "
                "BLQ concentrations were omitted and were not plotted as zero."
            ),
            "fits": [item.fit.model_dump(mode="json") for item in fits],
            "diagnostics": [item.diagnostic.model_dump(mode="json") for item in fits],
        },
    )


def _write_typed_parquet(
    path: Path,
    columns: tuple[str, ...],
    floats: set[str],
    rows: list[dict[str, object]],
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    fields = [pa.field(name, pa.float64() if name in floats else pa.string()) for name in columns]
    stored = [{column: row.get(column) for column in columns} for row in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(stored, schema=pa.schema(fields)), path)


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
