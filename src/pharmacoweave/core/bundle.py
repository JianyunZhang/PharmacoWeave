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

from pharmacoweave.core.entities import AssayQuantitation, ObservationType
from pharmacoweave.core.evidence import EvidenceExchange, EvidenceRecord, ReviewStatus
from pharmacoweave.core.modeling import NcaParameterSet
from pharmacoweave.core.provenance import (
    BundleFile,
    OmittedArtifact,
    Producer,
    ResultBundleManifest,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import CONC_AMOUNT, CONC_MASS, convert, parse_unit
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import IngestError, UnitError, ValidationFailed
from pharmacoweave.ingestion.tabular import (
    analyte_row,
    assay_row,
    compound_row,
    dose_row,
    identity_row,
    observation_row,
    sample_row,
    subject_row,
    urine_collection_row,
    write_named_table,
)
from pharmacoweave.reporting.figures import (
    EffectCurveOverlay,
    EffectPlotPoint,
    PredictedPoint,
    SemilogOverlay,
    concentration_effect_svg,
    concentration_time_semilog_svg,
    concentration_time_svg,
    goodness_of_fit_svg,
    residual_diagnostics_svg,
    weighted_residual_svg,
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
        role="derived_observations",
        reason="analysis/derived_observations.parquet is absent because no analysis or BLQ policy was run.",
    ),
    OmittedArtifact(
        role="nca_parameters",
        reason="analysis/nca_parameters.parquet is absent because NCA was not run.",
    ),
    OmittedArtifact(
        role="urine_excretion",
        reason="analysis/urine_excretion.parquet is absent because urine excretion NCA was not run.",
    ),
    OmittedArtifact(
        role="model_parameters",
        reason="analysis/model_parameters.parquet is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="predictions",
        reason="analysis/predictions.parquet is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="residuals",
        reason="analysis/residuals.parquet is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="emax_parameters",
        reason="analysis/emax_parameters.parquet is absent because no PD Emax model was fit.",
    ),
    OmittedArtifact(
        role="emax_predictions",
        reason="analysis/emax_predictions.parquet is absent because no PD Emax model was fit.",
    ),
    OmittedArtifact(
        role="diagnostics",
        reason="analysis/diagnostics.json is absent because no model was fit.",
    ),
    OmittedArtifact(
        role="concentration_time_figure",
        reason="figures/concentration_time.svg is absent because no PK analysis figure was drawn.",
    ),
    OmittedArtifact(
        role="concentration_time_semilog_figure",
        reason="figures/concentration_time_semilog.svg is absent because no PK analysis figure was drawn.",
    ),
    OmittedArtifact(
        role="goodness_of_fit_figure",
        reason="figures/goodness_of_fit.svg is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="residual_diagnostics_figure",
        reason="figures/residual_diagnostics.svg is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="weighted_residuals_figure",
        reason="figures/weighted_residuals.svg is absent because no PK model was fit.",
    ),
    OmittedArtifact(
        role="concentration_effect_figure",
        reason="figures/concentration_effect.svg is absent because no PD Emax model was fit.",
    ),
    OmittedArtifact(
        role="ivive_predictions",
        reason="analysis/ivive_predictions.parquet is absent because no exploratory IVIVE was run.",
    ),
    OmittedArtifact(
        role="prediction_audit",
        reason="analysis/prediction_audit.json is absent because no prior-vs-estimate audit was run.",
    ),
    OmittedArtifact(
        role="exposure_simulation",
        reason="analysis/exposure_simulation.json is absent because no exploratory exposure simulation was run.",
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
    emax_fits: Sequence[Any] | None = None,
    adme_records: Sequence[Any] | None = None,
    tox_endpoints: Sequence[Any] | None = None,
    blq_policy: Any | None = None,
    ivive_predictions: Sequence[Any] | None = None,
    prediction_audit: Any | None = None,
    exposure_simulation: Any | None = None,
    urine_results: Sequence[Any] | None = None,
) -> ResultBundleManifest:
    """Validate ``study`` and write a new bundle directory.

    Pass ``nca_parameters`` from ``run_nca`` to store NCA. Pass IV bolus
    ``fits`` to store PK parameters, residuals, and diagnostic figures. Pass
    ``emax_fits`` to store direct/sigmoid Emax curves, and ``adme_records`` or
    ``tox_endpoints`` to persist structured DMPK and toxicology records and
    screened ``PharmacologyEvidence``. Pass ``ivive_predictions`` (exploratory,
    predicted grade) and ``prediction_audit`` (a ``PredictionAuditReport``) to
    persist the IVIVE table and the prior-vs-estimate fold-error audit. Pass
    ``exposure_simulation`` (an ``ExposureSimulationResult``) to persist an
    exploratory single-dose curve. That curve stays predicted and is not a dose.
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
            emax_fits=emax_fits,
            adme_records=adme_records,
            tox_endpoints=tox_endpoints,
            blq_policy=blq_policy,
            ivive_predictions=ivive_predictions,
            prediction_audit=prediction_audit,
            exposure_simulation=exposure_simulation,
            urine_results=urine_results,
        )
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _write_bundle(
    study: PharmacologyStudy,
    destination: Path,
    *,
    run_id: str,
    report,
    nca_parameters: list[NcaParameterSet] | None,
    fits: Sequence[Any] | None,
    emax_fits: Sequence[Any] | None,
    adme_records: Sequence[Any] | None,
    tox_endpoints: Sequence[Any] | None,
    blq_policy: Any | None,
    ivive_predictions: Sequence[Any] | None = None,
    prediction_audit: Any | None = None,
    exposure_simulation: Any | None = None,
    urine_results: Sequence[Any] | None = None,
) -> ResultBundleManifest:
    created_at = datetime.now(UTC)
    counts = _counts(study, report, adme_records=adme_records, tox_endpoints=tox_endpoints)
    if ivive_predictions is not None:
        counts["n_ivive_predictions"] = len(ivive_predictions)
    if prediction_audit is not None:
        counts["n_prediction_audit_entries"] = len(prediction_audit.entries)
        counts["n_prediction_audit_skipped"] = len(prediction_audit.skipped)
    if exposure_simulation is not None:
        counts["n_simulated_times"] = len(exposure_simulation.curve)
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
    if study.urine_collections:
        write_named_table(
            destination / "experimental/urine_collections.parquet",
            "urine_collections",
            [urine_collection_row(item) for item in study.urine_collections],
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

    paths = dict(_PATHS)
    if study.urine_collections:
        paths["urine_collections"] = "experimental/urine_collections.parquet"
    if adme_records is not None:
        adme_path = "experimental/in_vitro_adme.parquet"
        paths["in_vitro_adme"] = adme_path
        _write_adme_parquet(destination / adme_path, adme_records)
    if tox_endpoints is not None:
        tox_path = "experimental/tox_endpoints.parquet"
        paths["tox_endpoints"] = tox_path
        _write_tox_parquet(destination / tox_path, tox_endpoints)

    _write_json(
        destination / _PATHS["source_files"],
        {"files": [item.model_dump(mode="json") for item in study.source_files]},
    )
    analysis_executed = any(
        item is not None
        for item in (
            nca_parameters,
            fits,
            emax_fits,
            blq_policy,
            ivive_predictions,
            prediction_audit,
            exposure_simulation,
            urine_results,
        )
    )
    notes: list[str] = []
    if nca_parameters is not None:
        notes.append("Single-dose NCA was run.")
    if fits is not None:
        notes.append("Compartmental PK fit was run.")
    if emax_fits is not None:
        notes.append("Direct/sigmoid Emax exposure-response fitting was run.")
    if ivive_predictions is not None:
        notes.append("Exploratory explicit-assumption well-stirred IVIVE was run; outputs stay predicted.")
    if prediction_audit is not None:
        notes.append("Prior predictions were audited against in vivo estimates with fold-error metrics.")
    if urine_results is not None:
        notes.append(
            "Urine excretion NCA was run. Amounts use recorded void volumes. "
            "Renal clearance uses a caller-supplied plasma AUC and is not plasma clearance."
        )
    if exposure_simulation is not None:
        notes.append(
            "Exploratory single-dose exposure simulation was run from declared parameters. "
            "The curve stays predicted and is not a dosing recommendation."
        )
    if blq_policy is not None and not notes:
        notes.append("Derived observations were materialized under an explicit BLQ policy.")

    estimation_methods: set[str] = set()
    if fits is not None:
        estimation_methods.update(item.fit.estimation_method for item in fits)
    if emax_fits is not None:
        estimation_methods.update(item.fit.estimation_method for item in emax_fits)

    _write_json(
        destination / _PATHS["analysis_parameters"],
        {
            "analysis_executed": analysis_executed,
            "parameters": {} if blq_policy is None else {"blq_policy": blq_policy.model_dump(mode="json")},
            "note": " ".join(notes) if notes else _ANALYSIS_NOTE,
            "terminal_phase_rule": None if nca_parameters is None else "pharmacoweave.lambda_z.adjusted_r2.v1",
            "estimation_methods": sorted(estimation_methods),
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

    exchange = _build_evidence_exchange(
        study,
        run_id=run_id,
        nca_parameters=nca_parameters,
        fits=fits,
        emax_fits=emax_fits,
        adme_records=adme_records,
        tox_endpoints=tox_endpoints,
        ivive_predictions=ivive_predictions,
    )
    counts["n_pharmacology_evidence"] = len(exchange.records)
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

    omitted = list(_OMITTED)
    if analysis_executed:
        from pharmacoweave.bioanalysis.qc import derive_study_observations

        derived_rows = derive_study_observations(study, policy=blq_policy)
        derived_path = "analysis/derived_observations.parquet"
        paths["derived_observations"] = derived_path
        _write_derived_observations_parquet(destination / derived_path, derived_rows)
        omitted = [item for item in omitted if item.role != "derived_observations"]

    semilog_text: str | None = None
    if nca_parameters is not None or fits is not None:
        from pharmacoweave.core.profiles import single_dose_profiles

        profiles = single_dose_profiles(study)
        semilog_overlays = _semilog_overlays(fits)
        term_counts = _terminal_point_counts(nca_parameters)
        semilog_text = concentration_time_semilog_svg(
            profiles,
            lloq_mg_l=_study_lloq_mg_l(study),
            overlays=semilog_overlays,
            terminal_point_counts=term_counts,
        )
        semilog_path = "figures/concentration_time_semilog.svg"
        paths["concentration_time_semilog_figure"] = semilog_path
        semilog_file = destination / semilog_path
        semilog_file.parent.mkdir(parents=True, exist_ok=True)
        semilog_file.write_text(semilog_text, encoding="utf-8")
        omitted = [item for item in omitted if item.role != "concentration_time_semilog_figure"]

        if nca_parameters is not None:
            nca_path = "analysis/nca_parameters.parquet"
            figure_path = "figures/concentration_time.svg"
            paths["nca_parameters"] = nca_path
            paths["concentration_time_figure"] = figure_path
            _write_nca_parquet(destination / nca_path, nca_parameters)
            figure = destination / figure_path
            figure.parent.mkdir(parents=True, exist_ok=True)
            figure.write_text(concentration_time_svg(profiles), encoding="utf-8")
            omitted = [item for item in omitted if item.role not in {"nca_parameters", "concentration_time_figure"}]

    if urine_results is not None:
        urine_path = "analysis/urine_excretion.parquet"
        paths["urine_excretion"] = urine_path
        _write_urine_parquet(destination / urine_path, urine_results)
        omitted = [item for item in omitted if item.role != "urine_excretion"]

    fit_points: list[PredictedPoint] = []
    gof_text: str | None = None
    residual_text: str | None = None
    wres_text: str | None = None
    if fits is not None:
        fit_points = _predicted_points(fits)
        fit_paths = {
            "model_parameters": "analysis/model_parameters.parquet",
            "predictions": "analysis/predictions.parquet",
            "residuals": "analysis/residuals.parquet",
        }
        paths.update(fit_paths)
        _write_fit_tables(destination, fits)
        omitted = [item for item in omitted if item.role not in fit_paths]
        if fit_points:
            gof_path = "figures/goodness_of_fit.svg"
            residual_path = "figures/residual_diagnostics.svg"
            wres_path = "figures/weighted_residuals.svg"
            paths["goodness_of_fit_figure"] = gof_path
            paths["residual_diagnostics_figure"] = residual_path
            paths["weighted_residuals_figure"] = wres_path
            gof_text = goodness_of_fit_svg(fit_points)
            residual_text = residual_diagnostics_svg(fit_points)
            wres_text = weighted_residual_svg(fit_points)
            gof_file = destination / gof_path
            gof_file.parent.mkdir(parents=True, exist_ok=True)
            gof_file.write_text(gof_text, encoding="utf-8")
            (destination / residual_path).write_text(residual_text, encoding="utf-8")
            (destination / wres_path).write_text(wres_text, encoding="utf-8")
            omitted = [
                item
                for item in omitted
                if item.role
                not in {"goodness_of_fit_figure", "residual_diagnostics_figure", "weighted_residuals_figure"}
            ]
        else:
            omitted = [_figure_omitted(item) for item in omitted]

    effect_svg_text: str | None = None
    if emax_fits is not None:
        emax_paths = {
            "emax_parameters": "analysis/emax_parameters.parquet",
            "emax_predictions": "analysis/emax_predictions.parquet",
        }
        paths.update(emax_paths)
        _write_emax_tables(destination, emax_fits)
        omitted = [item for item in omitted if item.role not in emax_paths]
        effect_svg_text = _render_emax_svg(emax_fits)
        if effect_svg_text is not None:
            eff_fig_path = "figures/concentration_effect.svg"
            paths["concentration_effect_figure"] = eff_fig_path
            eff_file = destination / eff_fig_path
            eff_file.parent.mkdir(parents=True, exist_ok=True)
            eff_file.write_text(effect_svg_text, encoding="utf-8")
            omitted = [item for item in omitted if item.role != "concentration_effect_figure"]

    if fits is not None or emax_fits is not None:
        diag_path = "analysis/diagnostics.json"
        paths["diagnostics"] = diag_path
        _write_diagnostics_json(destination / diag_path, fits=fits, emax_fits=emax_fits)
        omitted = [item for item in omitted if item.role != "diagnostics"]

    if ivive_predictions is not None:
        ivive_path = "analysis/ivive_predictions.parquet"
        paths["ivive_predictions"] = ivive_path
        _write_ivive_parquet(destination / ivive_path, ivive_predictions)
        omitted = [item for item in omitted if item.role != "ivive_predictions"]
    if prediction_audit is not None:
        audit_path = "analysis/prediction_audit.json"
        paths["prediction_audit"] = audit_path
        _write_json(destination / audit_path, json.loads(prediction_audit.model_dump_json()))
        omitted = [item for item in omitted if item.role != "prediction_audit"]
    if exposure_simulation is not None:
        simulation_path = "analysis/exposure_simulation.json"
        paths["exposure_simulation"] = simulation_path
        _write_json(destination / simulation_path, json.loads(exposure_simulation.model_dump_json()))
        omitted = [item for item in omitted if item.role != "exposure_simulation"]

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
        emax_rows=_emax_rows(emax_fits),
        semilog_svg=semilog_text,
        goodness_of_fit_svg=gof_text,
        residual_diagnostics_svg=residual_text,
        weighted_residuals_svg=wres_text,
        concentration_effect_svg=effect_svg_text,
        ivive_rows=_ivive_rows(ivive_predictions),
        audit_summary_rows=_audit_summary_rows(prediction_audit),
        audit_entry_rows=_audit_entry_rows(prediction_audit),
        audit_skipped_rows=_audit_skipped_rows(prediction_audit),
        simulation_rows=_simulation_rows(exposure_simulation),
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


def _build_evidence_exchange(
    study: PharmacologyStudy,
    *,
    run_id: str,
    nca_parameters: list[NcaParameterSet] | None,
    fits: Sequence[Any] | None,
    emax_fits: Sequence[Any] | None,
    adme_records: Sequence[Any] | None,
    tox_endpoints: Sequence[Any] | None,
    ivive_predictions: Sequence[Any] | None = None,
    urine_results: Sequence[Any] | None = None,
) -> EvidenceExchange:
    has_run_sources = any(
        item is not None
        for item in (nca_parameters, fits, emax_fits, adme_records, tox_endpoints, ivive_predictions, urine_results)
    )
    if not has_run_sources:
        return EvidenceExchange(records=list(study.external_evidence))

    from pharmacoweave.integrations.pharmacology_evidence import collect_run_evidence

    generated = collect_run_evidence(
        run_id=run_id,
        nca_parameters=nca_parameters,
        pk_fits=fits,
        emax_fits=emax_fits,
        adme_records=adme_records,
        tox_endpoints=tox_endpoints,
        ivive_predictions=ivive_predictions,
        urine_results=urine_results,
    )
    merged: list[EvidenceRecord] = list(study.external_evidence)
    seen = {record.evidence_id for record in merged}
    for record in generated:
        if record.evidence_id not in seen:
            merged.append(record)
            seen.add(record.evidence_id)
    return EvidenceExchange(records=merged)


def _counts(
    study: PharmacologyStudy,
    report,
    *,
    adme_records: Sequence[Any] | None = None,
    tox_endpoints: Sequence[Any] | None = None,
) -> dict[str, int]:
    observations = study.observations
    counts = {
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
    if adme_records is not None:
        counts["n_in_vitro_adme"] = len(adme_records)
    if tox_endpoints is not None:
        counts["n_tox_endpoints"] = len(tox_endpoints)
    return counts


def _study_lloq_mg_l(study: PharmacologyStudy) -> float | None:
    analytes = {item.analyte_id: item for item in study.analytes}
    compounds = {item.compound_id: item for item in study.compounds}
    lloqs: list[float] = []
    for assay in study.assays:
        if assay.quantitation is not AssayQuantitation.CONCENTRATION:
            continue
        if assay.lloq is None or assay.limit_unit is None:
            continue
        try:
            dims = parse_unit(assay.limit_unit).dimensions
            if dims == CONC_MASS:
                lloqs.append(convert(assay.lloq, assay.limit_unit, "mg/L"))
            elif dims == CONC_AMOUNT:
                analyte = analytes.get(assay.analyte_id) if assay.analyte_id is not None else None
                compound = (
                    compounds.get(analyte.compound_id)
                    if analyte is not None and analyte.compound_id is not None
                    else None
                )
                if compound is not None and compound.molecular_weight is not None:
                    lloqs.append(compound.convert_quantity(assay.lloq, assay.limit_unit, "mg/L"))
        except UnitError:
            continue
    return max(lloqs) if lloqs else None


def _terminal_point_counts(nca_parameters: list[NcaParameterSet] | None) -> dict[str, int]:
    if not nca_parameters:
        return {}
    counts: dict[str, int] = {}
    for row in nca_parameters:
        if row.terminal_phase_qualified and row.terminal_phase_n_points:
            counts[row.subject_id] = int(row.terminal_phase_n_points)
    return counts


def _semilog_overlays(fits: Sequence[Any] | None) -> list[SemilogOverlay]:
    if not fits:
        return []
    overlays: list[SemilogOverlay] = []
    for item in fits:
        predicted = tuple(item.predicted_mg_l)
        if predicted:
            overlays.append(
                SemilogOverlay(
                    series_id=item.subject_id,
                    times_h=tuple(item.times_h),
                    predicted_mg_l=predicted,
                )
            )
    return overlays


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
        "regimen",
        "n_doses",
        "dosing_interval_h",
        "auc_tau",
        "auc_tau_unit",
        "cmin",
        "cmin_unit",
        "cavg",
        "cavg_unit",
        "fluctuation",
        "accumulation_ratio",
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
        "dosing_interval_h",
        "auc_tau",
        "cmin",
        "cavg",
        "fluctuation",
        "accumulation_ratio",
    }
    ints = {"terminal_phase_n_points", "n_doses"}
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


def _write_derived_observations_parquet(path: Path, rows: Sequence[Any]) -> None:
    columns = (
        "observation_id",
        "sample_id",
        "subject_id",
        "analyte_id",
        "time_since_dose",
        "time_unit",
        "is_predose",
        "stored_status",
        "stored_numeric_value",
        "stored_unit",
        "included_in_analysis",
        "analysis_value",
        "analysis_unit",
        "policy_rule",
        "caution",
    )
    dicts = [item.model_dump(mode="json") for item in rows]
    _write_typed_parquet(
        path,
        columns,
        floats={"time_since_dose", "stored_numeric_value", "analysis_value"},
        bools={"is_predose", "included_in_analysis"},
        rows=dicts,
    )


def _write_adme_parquet(path: Path, records: Sequence[Any]) -> None:
    columns = (
        "record_id",
        "compound_id",
        "assay_kind",
        "system",
        "status",
        "value",
        "unit",
        "replicate_id",
        "qc_reason",
        "evidence_grade",
        "research_use_only",
        "limitations",
    )
    rows: list[dict[str, object]] = []
    for item in records:
        payload = item.model_dump(mode="json")
        payload["limitations"] = " | ".join(item.limitations)
        rows.append(payload)
    _write_typed_parquet(
        path,
        columns,
        floats={"value"},
        bools={"research_use_only"},
        rows=rows,
    )


def _write_tox_parquet(path: Path, endpoints: Sequence[Any]) -> None:
    columns = (
        "endpoint_id",
        "compound_id",
        "kind",
        "finding",
        "severity",
        "subject_id",
        "dose_event_id",
        "exposure_ref",
        "numeric_value",
        "unit",
        "incidence_affected",
        "incidence_examined",
        "evidence_grade",
        "research_use_only",
        "limitations",
    )
    rows: list[dict[str, object]] = []
    for item in endpoints:
        payload = item.model_dump(mode="json")
        payload["limitations"] = " | ".join(item.limitations)
        rows.append(payload)
    _write_typed_parquet(
        path,
        columns,
        floats={"numeric_value"},
        ints={"incidence_affected", "incidence_examined"},
        bools={"research_use_only"},
        rows=rows,
    )


def _write_ivive_parquet(path: Path, records: Sequence[Any]) -> None:
    columns = (
        "prediction_id",
        "compound_id",
        "species",
        "physiology_profile_id",
        "source_record_id",
        "assay_kind",
        "clint_app_value",
        "clint_app_unit",
        "fu_inc",
        "fu_plasma",
        "blood_to_plasma_ratio",
        "fu_blood",
        "clint_u_vitro_value",
        "scaling_factor_value",
        "scaling_factor_unit",
        "liver_weight_g",
        "body_weight_kg",
        "hepatic_blood_flow_l_h",
        "clint_u_hepatic_l_h",
        "hepatic_extraction_ratio",
        "hepatic_blood_clearance_l_h",
        "hepatic_plasma_clearance_l_h",
        "hepatic_plasma_clearance_l_h_kg",
        "model",
        "clearance_unit",
        "evidence_grade",
        "research_use_only",
        "assumptions",
        "limitations",
        "derived_from",
    )
    rows: list[dict[str, object]] = []
    for item in records:
        payload = item.model_dump(mode="json")
        payload["assumptions"] = " | ".join(item.assumptions)
        payload["limitations"] = " | ".join(item.limitations)
        payload["derived_from"] = " | ".join(item.derived_from)
        rows.append(payload)
    _write_typed_parquet(
        path,
        columns,
        floats={
            "clint_app_value",
            "fu_inc",
            "fu_plasma",
            "blood_to_plasma_ratio",
            "fu_blood",
            "clint_u_vitro_value",
            "scaling_factor_value",
            "liver_weight_g",
            "body_weight_kg",
            "hepatic_blood_flow_l_h",
            "clint_u_hepatic_l_h",
            "hepatic_extraction_ratio",
            "hepatic_blood_clearance_l_h",
            "hepatic_plasma_clearance_l_h",
            "hepatic_plasma_clearance_l_h_kg",
        },
        bools={"research_use_only"},
        rows=rows,
    )


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
    if item.role == "weighted_residuals_figure":
        return OmittedArtifact(
            role=item.role,
            reason="figures/weighted_residuals.svg is absent because the fit produced no positive predictions.",
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
        by_param = {estimate.parameter: estimate for estimate in item.estimates}
        cl_est = by_param.get("CL") or by_param.get("CL/F")
        v_est = by_param.get("V") or by_param.get("V/F") or by_param.get("V1")
        summary_parts = [
            f"{estimate.parameter} = {_fmt_param_with_unit(estimate)}"
            for estimate in item.estimates
        ]
        rows.append(
            {
                "subject_id": item.subject_id,
                "model_id": item.fit.model_id,
                "method": item.fit.estimation_method,
                "status": item.diagnostic.status.value,
                "parameters_summary": "; ".join(summary_parts),
                "clearance": _fmt_fit(None if cl_est is None else cl_est.value),
                "clearance_ci": _fmt_ci(cl_est),
                "volume": _fmt_fit(None if v_est is None else v_est.value),
                "volume_ci": _fmt_ci(v_est),
            }
        )
    return rows


def _emax_rows(emax_fits: Sequence[Any] | None) -> list[dict[str, str]]:
    if not emax_fits:
        return []
    rows: list[dict[str, str]] = []
    for item in emax_fits:
        by_name = {est.name: est for est in item.estimates}
        rows.append(
            {
                "fit_id": item.fit.fit_id,
                "method": item.fit.estimation_method,
                "status": item.diagnostic.status.value,
                "e0": _fmt_param_with_unit(by_name.get("E0")),
                "emax": _fmt_param_with_unit(by_name.get("Emax")),
                "ec50": _fmt_param_with_unit(by_name.get("EC50")),
                "hill": _fmt_fit(None if "hill" not in by_name else by_name["hill"].value),
            }
        )
    return rows


def _simulation_rows(result: Any | None) -> list[dict[str, str]]:
    if result is None:
        return []
    threshold = "not requested"
    if result.time_above_threshold_h is not None:
        threshold = f"{result.time_above_threshold_h:.6g} h"
    return [
        {
            "simulation_id": result.simulation_id,
            "model_id": result.model_id,
            "cmax": f"{result.cmax:.6g} {result.concentration_unit}",
            "tmax_h": f"{result.tmax_h:.6g}",
            "auc_inf": f"{result.auc_inf:.6g} {result.auc_unit}",
            "half_life_h": f"{result.half_life_h:.6g}",
            "time_above": threshold,
            "evidence_type": result.evidence_type,
            "dosing_recommendation": "false",
        }
    ]


def _ivive_rows(predictions: Sequence[Any] | None) -> list[dict[str, str]]:
    if not predictions:
        return []
    rows: list[dict[str, str]] = []
    for item in predictions:
        rows.append(
            {
                "prediction_id": item.prediction_id,
                "compound_id": item.compound_id,
                "species": item.species,
                "assay_kind": item.assay_kind,
                "clint_app": f"{item.clint_app_value:.6g} {item.clint_app_unit}",
                "fu_plasma": f"{item.fu_plasma:.4g}",
                "fu_inc": f"{item.fu_inc:.4g}",
                "blood_to_plasma_ratio": f"{item.blood_to_plasma_ratio:.4g}",
                "cl_plasma_l_h": f"{item.hepatic_plasma_clearance_l_h:.6g}",
                "cl_plasma_l_h_kg": f"{item.hepatic_plasma_clearance_l_h_kg:.6g}",
                "extraction_ratio": f"{item.hepatic_extraction_ratio:.4g}",
            }
        )
    return rows


def _audit_summary_rows(report: Any | None) -> list[dict[str, str]]:
    if report is None or not report.summaries:
        return []
    return [
        {
            "measurement_type": s.measurement_type,
            "n_comparisons": str(s.n_comparisons),
            "afe": f"{s.afe:.4g}",
            "aafe": f"{s.aafe:.4g}",
            "fraction_within_2_fold": f"{s.fraction_within_2_fold * 100:.1f}%",
            "fraction_within_3_fold": f"{s.fraction_within_3_fold * 100:.1f}%",
            "n_outlier_2_fold": str(s.n_outlier_2_fold),
        }
        for s in report.summaries
    ]


def _audit_entry_rows(report: Any | None) -> list[dict[str, str]]:
    if report is None or not report.entries:
        return []
    return [
        {
            "prior_id": e.prior_evidence_id,
            "estimate_id": e.estimate_evidence_id,
            "measurement_type": e.measurement_type,
            "predicted": f"{e.predicted_value:.6g} {e.unit}",
            "estimated": f"{e.estimated_value:.6g} {e.unit}",
            "ratio": f"{e.ratio:.4g}",
            "fold_error": f"{e.fold_error:.4g}",
            "bias": e.bias_direction.value,
        }
        for e in report.entries
    ]


def _audit_skipped_rows(report: Any | None) -> list[dict[str, str]]:
    if report is None or not report.skipped:
        return []
    return [
        {
            "prior_id": s.prior_evidence_id,
            "reason": s.reason,
        }
        for s in report.skipped
    ]


def _fmt_param_with_unit(estimate: Any | None) -> str:
    if estimate is None:
        return ""
    base = f"{estimate.value:.6g} {estimate.unit}"
    ci = _fmt_ci(estimate)
    return f"{base} {ci}" if ci else base


def _fmt_ci(estimate: Any | None) -> str:
    if estimate is None or estimate.ci_lower is None or estimate.ci_upper is None:
        return ""
    return f"[{estimate.ci_lower:.4g}, {estimate.ci_upper:.4g}]"


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
                    "standard_error": estimate.standard_error,
                    "ci_lower": estimate.ci_lower,
                    "ci_upper": estimate.ci_upper,
                    "ci_level": estimate.ci_level,
                    "bound_hit": estimate.bound_hit,
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
            rel_residual = residual / pred
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
                    "relative_residual": rel_residual,
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
            "standard_error",
            "ci_lower",
            "ci_upper",
            "ci_level",
            "bound_hit",
        ),
        floats={"value", "standard_error", "ci_lower", "ci_upper", "ci_level"},
        bools={"bound_hit"},
        rows=parameter_rows,
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
        floats={"time_h", "observed_mg_l", "predicted_mg_l"},
        rows=prediction_rows,
    )
    _write_typed_parquet(
        destination / "analysis" / "residuals.parquet",
        ("fit_id", "observation_id", "time_h", "residual_mg_l", "relative_residual", "residual_unit"),
        floats={"time_h", "residual_mg_l", "relative_residual"},
        rows=residual_rows,
    )


def _write_emax_tables(destination: Path, emax_fits: Sequence[Any]) -> None:
    param_rows: list[dict[str, object]] = []
    pred_rows: list[dict[str, object]] = []
    for item in emax_fits:
        for est in item.estimates:
            param_rows.append(
                {
                    "fit_id": est.fit_id,
                    "name": est.name,
                    "value": est.value,
                    "unit": est.unit,
                    "fixed": est.fixed,
                    "estimation_method": est.estimation_method,
                    "standard_error": est.standard_error,
                    "ci_lower": est.ci_lower,
                    "ci_upper": est.ci_upper,
                    "ci_level": est.ci_level,
                }
            )
        pair_by_effect = {pair.effect_id: pair for pair in item.pairs}
        for pred in item.predicted:
            pair = pair_by_effect.get(pred.effect_id)
            pred_rows.append(
                {
                    "fit_id": item.fit.fit_id,
                    "effect_id": pred.effect_id,
                    "exposure_ref": pred.exposure_ref,
                    "concentration": None if pair is None else pair.concentration,
                    "concentration_unit": None if pair is None else pair.concentration_unit,
                    "observed_effect": None if pair is None else pair.effect,
                    "predicted_effect": pred.value,
                    "effect_unit": pred.unit,
                    "evidence_grade": pred.evidence_grade,
                }
            )
    _write_typed_parquet(
        destination / "analysis" / "emax_parameters.parquet",
        (
            "fit_id",
            "name",
            "value",
            "unit",
            "fixed",
            "estimation_method",
            "standard_error",
            "ci_lower",
            "ci_upper",
            "ci_level",
        ),
        floats={"value", "standard_error", "ci_lower", "ci_upper", "ci_level"},
        bools={"fixed"},
        rows=param_rows,
    )
    _write_typed_parquet(
        destination / "analysis" / "emax_predictions.parquet",
        (
            "fit_id",
            "effect_id",
            "exposure_ref",
            "concentration",
            "concentration_unit",
            "observed_effect",
            "predicted_effect",
            "effect_unit",
            "evidence_grade",
        ),
        floats={"concentration", "observed_effect", "predicted_effect"},
        rows=pred_rows,
    )


def _render_emax_svg(emax_fits: Sequence[Any]) -> str | None:
    points: list[EffectPlotPoint] = []
    overlays: list[EffectCurveOverlay] = []
    conc_unit = "mg/L"
    eff_unit = "effect"
    for item in emax_fits:
        if item.pairs:
            conc_unit = item.pairs[0].concentration_unit
            eff_unit = item.pairs[0].effect_unit
        for pair in item.pairs:
            points.append(
                EffectPlotPoint(
                    series_id=item.fit.fit_id,
                    concentration=pair.concentration,
                    effect=pair.effect,
                )
            )
        by_name = {est.name: est.value for est in item.estimates}
        if {"E0", "Emax", "EC50"} <= set(by_name):
            overlays.append(
                EffectCurveOverlay(
                    series_id=item.fit.fit_id,
                    e0=by_name["E0"],
                    emax=by_name["Emax"],
                    ec50=by_name["EC50"],
                    hill=by_name.get("hill", 1.0),
                )
            )
    if not points:
        return None
    return concentration_effect_svg(
        points,
        overlays=overlays,
        concentration_unit=conc_unit,
        effect_unit=eff_unit,
    )


def _write_diagnostics_json(
    path: Path,
    *,
    fits: Sequence[Any] | None,
    emax_fits: Sequence[Any] | None,
) -> None:
    all_items = [*(fits or ()), *(emax_fits or ())]
    _write_json(
        path,
        {
            "research_use_only": True,
            "note": (
                "A passing diagnostic states limitations. "
                "BLQ concentrations and censored effects were omitted and were not plotted as zero."
            ),
            "fits": [item.fit.model_dump(mode="json") for item in all_items],
            "diagnostics": [item.diagnostic.model_dump(mode="json") for item in all_items],
        },
    )


def _write_urine_parquet(path: Path, results: Sequence[Any]) -> None:
    columns = (
        "subject_id",
        "analyte_id",
        "dose_event_id",
        "matrix",
        "observation_id",
        "interval_start_h",
        "interval_end_h",
        "urine_concentration_mg_l",
        "urine_volume_l",
        "amount_excreted_mg",
        "omitted_reason",
        "total_amount_excreted_mg",
        "renal_clearance_l_h",
        "plasma_auc_source_id",
    )
    floats = {
        "interval_start_h",
        "interval_end_h",
        "urine_concentration_mg_l",
        "urine_volume_l",
        "amount_excreted_mg",
        "total_amount_excreted_mg",
        "renal_clearance_l_h",
    }
    rows: list[dict[str, object]] = []
    for result in results:
        for interval in result.intervals:
            rows.append(
                {
                    "subject_id": result.subject_id,
                    "analyte_id": result.analyte_id,
                    "dose_event_id": result.dose_event_id,
                    "matrix": result.matrix,
                    "observation_id": interval.observation_id,
                    "interval_start_h": interval.interval_start_h,
                    "interval_end_h": interval.interval_end_h,
                    "urine_concentration_mg_l": interval.urine_concentration_mg_l,
                    "urine_volume_l": interval.urine_volume_l,
                    "amount_excreted_mg": interval.amount_excreted_mg,
                    "omitted_reason": interval.omitted_reason,
                    "total_amount_excreted_mg": result.amount_excreted_mg,
                    "renal_clearance_l_h": result.renal_clearance_l_h,
                    "plasma_auc_source_id": result.plasma_auc_source_id,
                }
            )
    _write_typed_parquet(path, columns, floats, rows)


def _write_typed_parquet(
    path: Path,
    columns: tuple[str, ...],
    floats: set[str],
    rows: list[dict[str, object]],
    *,
    ints: set[str] | None = None,
    bools: set[str] | None = None,
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    int_cols = ints or set()
    bool_cols = bools or set()
    fields: list[pa.Field] = []
    for name in columns:
        if name in floats:
            dtype: pa.DataType = pa.float64()
        elif name in int_cols:
            dtype = pa.int64()
        elif name in bool_cols:
            dtype = pa.bool_()
        else:
            dtype = pa.string()
        fields.append(pa.field(name, dtype))
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
