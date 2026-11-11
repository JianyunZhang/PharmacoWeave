"""Command line for validation, NCA, PK/PD fitting, schema export, and capabilities."""

import argparse
import json
from pathlib import Path
from typing import Any, Literal

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.errors import AnalysisError, IngestError, PharmacoWeaveError, ValidationFailed
from pharmacoweave.version import __version__


def main(argv: list[str] | None = None) -> int:
    """Run the PharmacoWeave command line. Returns a process status code."""

    parser = argparse.ArgumentParser(prog="pharmacoweave")
    parser.add_argument("--version", action="version", version=f"pharmacoweave {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    validate = commands.add_parser("validate", help="Validate a study directory")
    validate.add_argument("study_dir", type=Path)
    validate.add_argument("--strict", action="store_true", help="Promote warnings to errors")
    validate.add_argument("--bundle", type=Path, help="Write a validation bundle to this new directory")
    validate.add_argument("--run-id", help="Run identifier used when --bundle is set")

    nca = commands.add_parser("nca", help="Run single-dose rich-serial NCA")
    nca.add_argument("study_dir", type=Path)
    nca.add_argument("--config", type=Path, required=True)
    nca.add_argument(
        "--with-fit",
        action="store_true",
        help="Also run the matching 1-compartment PK fit alongside NCA when exporting a bundle",
    )
    nca.add_argument(
        "--bundle",
        type=Path,
        help="Write a bundle with NCA tables and concentration-time figures",
    )
    nca.add_argument("--run-id")

    urine = commands.add_parser("urine-nca", help="Run urine excretion amount and renal clearance")
    urine.add_argument("study_dir", type=Path)
    urine.add_argument("--config", type=Path, required=True)
    urine.add_argument("--plasma-auc", type=Path, help="YAML file of caller-supplied plasma AUC values")
    urine.add_argument("--bundle", type=Path, help="Write a bundle that includes urine excretion results")
    urine.add_argument("--run-id")

    fit = commands.add_parser(
        "fit",
        help="Fit a compartmental PK model (1c IV bolus, 1c infusion, 1c extravascular, 2c or 3c IV bolus)",
    )
    fit.add_argument("study_dir", type=Path)
    fit.add_argument(
        "--model",
        choices=(
            "one_compartment_iv_bolus",
            "one_compartment_iv_infusion",
            "one_compartment_extravascular",
            "two_compartment_iv_bolus",
            "three_compartment_iv_bolus",
        ),
        default="one_compartment_iv_bolus",
    )
    fit.add_argument("--method", choices=("log_linear_ols", "weighted_nls"), default=None)
    fit.add_argument("--weight", choices=("uniform", "1/y", "1/y^2"), default="1/y^2")
    fit.add_argument("--bundle", type=Path, help="Write a bundle that includes the fit and residual figures")
    fit.add_argument("--run-id")

    emax = commands.add_parser("emax", help="Fit a direct or sigmoid Emax curve on explicit concentrations")
    emax.add_argument("input_json", type=Path, help="JSON file conforming to DirectEmaxRequest")
    emax.add_argument("--method", choices=("nonlinear_ls", "linearized"), default="nonlinear_ls")
    emax.add_argument("--estimate-hill", action="store_true", help="Estimate the Sigmoid Hill coefficient")
    emax.add_argument("--evidence-out", type=Path, help="Write screened PharmacologyEvidence JSON to this path")
    emax.add_argument("--study", type=Path, help="Study directory used when exporting a result bundle")
    emax.add_argument("--bundle", type=Path, help="Write a result bundle including Emax outputs (requires --study)")
    emax.add_argument("--run-id")

    ivive = commands.add_parser("ivive", help="Run exploratory well-stirred hepatic IVIVE")
    ivive.add_argument("input_json", type=Path, help="JSON file containing an IviveRequest or list of IviveRequests")
    ivive.add_argument("--evidence-out", type=Path, help="Write screened PharmacologyEvidence JSON to this path")
    ivive.add_argument("--study", type=Path, help="Study directory used when exporting a result bundle")
    ivive.add_argument("--bundle", type=Path, help="Write a result bundle including IVIVE outputs (requires --study)")
    ivive.add_argument("--run-id")

    audit = commands.add_parser("audit-predictions", help="Audit prior predictions against in vivo estimates")
    audit.add_argument("study_dir", type=Path, help="Study directory")
    audit.add_argument(
        "--priors", type=Path, help="JSON file with prior EvidenceRecords (defaults to study.external_evidence)"
    )
    audit.add_argument("--nca-config", type=Path, help="NCA config YAML to compute in vivo estimates")
    audit.add_argument(
        "--fit-model",
        choices=(
            "one_compartment_iv_bolus",
            "one_compartment_iv_infusion",
            "one_compartment_extravascular",
            "two_compartment_iv_bolus",
            "three_compartment_iv_bolus",
        ),
        default=None,
        help="PK model to compute in vivo estimates",
    )
    audit.add_argument("--estimates", type=Path, help="JSON file with in vivo estimated EvidenceRecords")
    audit.add_argument("--output", type=Path, help="Write PredictionAuditReport JSON to this path")
    audit.add_argument("--bundle", type=Path, help="Write a result bundle with prediction audit")
    audit.add_argument("--run-id")

    verify = commands.add_parser(
        "verify-bundle",
        help="Verify ResultBundle cryptographic hashes, structure, and JSON schemas",
    )
    verify.add_argument("bundle_dir", type=Path, help="Path to ResultBundle directory")
    verify.add_argument("--json", action="store_true", help="Output verification report as JSON")

    diff = commands.add_parser(
        "diff",
        help="Compare two ResultBundles (manifest, file tree, parameters, and tables)",
    )
    diff.add_argument("bundle_a", type=Path, help="Path to first ResultBundle directory")
    diff.add_argument("bundle_b", type=Path, help="Path to second ResultBundle directory")
    diff.add_argument("--json", action="store_true", help="Output bundle diff as JSON")
    diff.add_argument(
        "--rel-tol",
        type=float,
        default=1e-5,
        help="Relative numerical tolerance for parameter comparison",
    )

    ingest = commands.add_parser(
        "ingest",
        help="Ingest raw tabular lab data using an ingest_mapping.yaml specification",
    )
    ingest.add_argument("--mapping", type=Path, required=True, help="Path to ingest_mapping.yaml")
    ingest.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Target directory where canonical study files will be written",
    )
    ingest.add_argument(
        "--format",
        choices=("csv", "tsv", "parquet"),
        default="csv",
        help="File format for exported study tables (default: csv)",
    )
    ingest.add_argument("--strict", action="store_true", help="Promote validation warnings to errors")
    ingest.add_argument("--bundle", type=Path, help="Also export a ResultBundle to this directory")
    ingest.add_argument("--run-id", help="Run identifier if exporting a bundle")

    simulate = commands.add_parser(
        "simulate",
        help="Simulate a single-dose analytical exposure interval from declared parameters",
    )
    simulate.add_argument("input_json", type=Path, help="JSON file conforming to ExposureSimulationRequest")
    simulate.add_argument("--output", type=Path, help="Write the ExposureSimulationResult JSON to this path")
    simulate.add_argument("--study", type=Path, help="Study directory used when exporting a result bundle")
    simulate.add_argument(
        "--bundle",
        type=Path,
        help="Write a result bundle including the simulation (requires --study)",
    )
    simulate.add_argument("--run-id")

    engine = commands.add_parser(
        "engine-run",
        help="Run the MIT analytical reference engine out of process, or refuse a copyleft engine",
    )
    engine.add_argument("input_json", type=Path, help="JSON file conforming to ExposureSimulationRequest")
    engine.add_argument(
        "--engine",
        default="pharmacoweave-analytical",
        help="Engine name. nlmixr2 and pk-sim are catalogued and refused.",
    )
    engine.add_argument("--output", type=Path, help="Write the EngineResponse JSON to this path")

    slot = commands.add_parser(
        "engine-slot",
        help="Run one caller-supplied executable when its file, version, and license text are all present",
    )
    slot.add_argument("request_file", type=Path, help="Request text passed to the executable")
    slot.add_argument("--executable", type=Path, required=True)
    slot.add_argument("--engine-version", required=True)
    slot.add_argument("--license-file", type=Path, required=True)
    slot.add_argument("--output", type=Path, help="Write the CallerEngineRun JSON to this path")

    accepted = commands.add_parser(
        "accepted-record",
        help="Evaluate one already accepted record and optionally store it in a result bundle",
    )
    accepted.add_argument("request_json", type=Path, help="JSON request with a kind field")
    accepted.add_argument("--output", type=Path, help="Write the function result JSON to this path")
    accepted.add_argument("--study", type=Path, help="Study directory used when exporting a result bundle")
    accepted.add_argument("--bundle", type=Path, help="Write a result bundle including the record (requires --study)")
    accepted.add_argument("--run-id")

    commands.add_parser("capabilities", help="Show which modules can execute")
    export = commands.add_parser("export-schemas", help="Write JSON Schema files from the data contracts")
    export.add_argument("--directory", type=Path, default=Path("schemas"))

    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _validate(args.study_dir, strict=args.strict, bundle=args.bundle, run_id=args.run_id)
        if args.command == "urine-nca":
            return _urine_nca(
                args.study_dir,
                config=args.config,
                plasma_auc=args.plasma_auc,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "nca":
            return _nca(
                args.study_dir,
                config=args.config,
                with_fit=args.with_fit,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "fit":
            return _fit(
                args.study_dir,
                model=args.model,
                method=args.method,
                weight=args.weight,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "accepted-record":
            return _accepted_record(
                args.request_json,
                output=args.output,
                study_dir=args.study,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "emax":
            return _emax(
                args.input_json,
                method=args.method,
                estimate_hill_coefficient=args.estimate_hill,
                evidence_out=args.evidence_out,
                study_dir=args.study,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "ivive":
            return _ivive(
                args.input_json,
                evidence_out=args.evidence_out,
                study_dir=args.study,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "audit-predictions":
            return _audit_predictions(
                args.study_dir,
                priors_path=args.priors,
                nca_config=args.nca_config,
                fit_model=args.fit_model,
                estimates_path=args.estimates,
                output=args.output,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "verify-bundle":
            return _verify_bundle(args.bundle_dir, as_json=args.json)
        if args.command == "diff":
            return _diff(args.bundle_a, args.bundle_b, as_json=args.json, rel_tol=args.rel_tol)
        if args.command == "ingest":
            return _ingest(
                args.mapping,
                output_dir=args.output_dir,
                file_format=args.format,
                strict=args.strict,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "simulate":
            return _simulate(
                args.input_json,
                output=args.output,
                study_dir=args.study,
                bundle=args.bundle,
                run_id=args.run_id,
            )
        if args.command == "engine-run":
            return _engine_run(args.input_json, engine=args.engine, output=args.output)
        if args.command == "engine-slot":
            return _engine_slot(
                args.request_file,
                executable=args.executable,
                engine_version=args.engine_version,
                license_file=args.license_file,
                output=args.output,
            )
        if args.command == "capabilities":
            return _capabilities()
        if args.command == "export-schemas":
            return _export_schemas(args.directory)
    except (IngestError, ValidationFailed, PharmacoWeaveError) as exc:
        print(exc)
        return 2
    return 2


def _validate(study_dir: Path, *, strict: bool, bundle: Path | None, run_id: str | None) -> int:
    from pharmacoweave.core.study import PharmacologyStudy

    mode = "strict" if strict else "research"
    study = PharmacologyStudy.from_files(study_dir, validate_schema=False)
    report = study.check(mode=mode)
    print(report)
    if not report.ok:
        return 1
    if bundle is not None:
        selected_run = run_id or f"{study.metadata.study_id}-validation"
        study.export_bundle(bundle, run_id=selected_run, mode=mode)
        print(f"Wrote validation bundle to {bundle}")
    return 0


def _nca(
    study_dir: Path,
    *,
    config: Path,
    with_fit: bool,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.core.modeling import ParameterSemantic
    from pharmacoweave.core.study import PharmacologyStudy
    from pharmacoweave.pk_nca.runner import load_nca_config

    loaded_cfg = load_nca_config(config)
    study = PharmacologyStudy.from_files(study_dir)
    results = study.run_nca(loaded_cfg)
    if with_fit and loaded_cfg.dosing != "single_dose":
        raise AnalysisError("--with-fit applies to single-dose NCA only.")
    for row in results:
        cl_label = "CL/F" if row.clearance_semantic is ParameterSemantic.CL_OVER_F else "CL"
        vz_label = "Vz/F" if row.volume_semantic is ParameterSemantic.VZ_OVER_F else "Vz"
        print(
            f"{row.subject_id}\tregimen={row.regimen}\tAUC_inf={row.auc_inf}\tAUC_tau={row.auc_tau}\t"
            f"{cl_label}={row.clearance}\thalf_life={row.half_life}\t{vz_label}={row.volume}"
        )
    if bundle is not None:
        fits: list[Any] | None = None
        if with_fit:
            if loaded_cfg.design_route == "iv_bolus":
                from pharmacoweave.pk_models.fit_iv import fit_iv_bolus

                fits = list(fit_iv_bolus(study))
            elif loaded_cfg.design_route in {"oral", "subcutaneous", "intramuscular"}:
                from pharmacoweave.pk_models.fit_extravascular import fit_one_compartment_extravascular

                fits = list(fit_one_compartment_extravascular(study))
            else:
                raise AnalysisError(
                    f"--with-fit is not implemented for route {loaded_cfg.design_route}; run NCA without --with-fit."
                )
        selected_run = run_id or f"{study.metadata.study_id}-nca"
        study.export_bundle(bundle, run_id=selected_run, nca_parameters=results, fits=fits)
        print(f"Wrote NCA bundle to {bundle}")
    return 0


def _urine_nca(
    study_dir: Path,
    *,
    config: Path,
    plasma_auc: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.core.study import PharmacologyStudy
    from pharmacoweave.pk_nca.urine import load_plasma_auc_inputs, load_urine_nca_config, run_urine_excretion

    loaded = load_urine_nca_config(config)
    aucs = load_plasma_auc_inputs(plasma_auc) if plasma_auc is not None else []
    study = PharmacologyStudy.from_files(study_dir)
    results = run_urine_excretion(study, loaded, plasma_aucs=aucs)
    for row in results:
        print(
            f"{row.subject_id}\tamount_excreted_mg={row.amount_excreted_mg}\t"
            f"renal_clearance_l_h={row.renal_clearance_l_h}"
        )
    if bundle is not None:
        selected_run = run_id or f"{study.metadata.study_id}-urine-nca"
        study.export_bundle(bundle, run_id=selected_run, urine_results=results)
        print(f"Wrote urine excretion bundle to {bundle}")
    return 0


def _fit(
    study_dir: Path,
    *,
    model: str,
    method: str | None,
    weight: str,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.core.study import PharmacologyStudy

    study = PharmacologyStudy.from_files(study_dir)
    if model == "one_compartment_iv_bolus":
        from pharmacoweave.pk_models.fit_iv import fit_iv_bolus
        from pharmacoweave.pk_models.wnls import fit_iv_bolus_wnls

        resolved_method = method or "log_linear_ols"
        fits: list[Any] = (
            list(fit_iv_bolus(study))
            if resolved_method == "log_linear_ols"
            else list(fit_iv_bolus_wnls(study, weight=weight))
        )
    elif model == "one_compartment_iv_infusion":
        from pharmacoweave.pk_models.fit_infusion import fit_one_compartment_iv_infusion

        if method == "log_linear_ols":
            raise AnalysisError("one_compartment_iv_infusion requires --method weighted_nls.")
        fits = list(fit_one_compartment_iv_infusion(study, weight=weight))
    elif model == "one_compartment_extravascular":
        from pharmacoweave.pk_models.fit_extravascular import fit_one_compartment_extravascular

        if method == "log_linear_ols":
            raise AnalysisError("one_compartment_extravascular requires --method weighted_nls.")
        fits = list(fit_one_compartment_extravascular(study, weight=weight))
    elif model == "two_compartment_iv_bolus":
        from pharmacoweave.pk_models.fit_two_compartment import fit_two_compartment_iv_bolus

        if method == "log_linear_ols":
            raise AnalysisError("two_compartment_iv_bolus requires --method weighted_nls.")
        fits = list(fit_two_compartment_iv_bolus(study, weight=weight))
    elif model == "three_compartment_iv_bolus":
        from pharmacoweave.pk_models.fit_three_compartment import fit_three_compartment_iv_bolus

        if method == "log_linear_ols":
            raise AnalysisError("three_compartment_iv_bolus requires --method weighted_nls.")
        fits = list(fit_three_compartment_iv_bolus(study, weight=weight))
    else:
        raise AnalysisError(f"Unsupported fit model {model}.")

    for item in fits:
        params = "\t".join(f"{est.parameter}={est.value:.6g}" for est in item.estimates)
        suffix = f"\t{params}" if params else ""
        print(
            f"{item.subject_id}\tmodel={item.fit.model_id}\t"
            f"method={item.fit.estimation_method}\tstatus={item.diagnostic.status.value}{suffix}"
        )
    if bundle is not None:
        selected_run = run_id or f"{study.metadata.study_id}-fit"
        study.export_bundle(bundle, run_id=selected_run, fits=fits)
        print(f"Wrote fit bundle to {bundle}")
    return 0


def _emax(
    input_json: Path,
    *,
    method: Literal["nonlinear_ls", "linearized"],
    estimate_hill_coefficient: bool,
    evidence_out: Path | None,
    study_dir: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.integrations.pharmacology_evidence import (
        evidence_from_emax,
        exchange_from_records,
        write_exchange,
    )
    from pharmacoweave.pkpd import DirectEmaxRequest, fit_linked_direct_emax

    if bundle is not None and study_dir is None:
        raise IngestError("--bundle requires --study so study provenance and validation are recorded.")
    if not input_json.is_file():
        raise IngestError(f"{input_json} does not exist.")
    try:
        payload = json.loads(input_json.read_text(encoding="utf-8"))
        request = DirectEmaxRequest.model_validate(payload)
    except Exception as exc:
        raise IngestError(f"{input_json} is not a valid DirectEmaxRequest: {exc}") from exc

    _linked, fit = fit_linked_direct_emax(
        request.effects,
        request.exposures,
        e0=request.e0,
        e0_unit=request.e0_unit,
        fit_id=request.fit_id,
        method=method,
        estimate_hill_coefficient=estimate_hill_coefficient,
    )
    by_name = {est.name: est.value for est in fit.estimates}
    emax_str = "" if "Emax" not in by_name else f"{by_name['Emax']:.6g}"
    ec50_str = "" if "EC50" not in by_name else f"{by_name['EC50']:.6g}"
    hill_str = "" if "hill" not in by_name else f"{by_name['hill']:.6g}"
    print(
        f"{fit.fit.fit_id}\tmethod={fit.fit.estimation_method}\t"
        f"status={fit.diagnostic.status.value}\tEmax={emax_str}\tEC50={ec50_str}\thill={hill_str}"
    )

    selected_run = run_id or f"{request.fit_id}-emax"
    if evidence_out is not None:
        exchange = exchange_from_records(evidence_from_emax(fit, run_id=selected_run))
        evidence_out.parent.mkdir(parents=True, exist_ok=True)
        write_exchange(exchange, evidence_out)
        print(f"Wrote PharmacologyEvidence to {evidence_out}")

    if bundle is not None and study_dir is not None:
        from pharmacoweave.core.study import PharmacologyStudy

        study = PharmacologyStudy.from_files(study_dir)
        study.export_bundle(bundle, run_id=selected_run, emax_fits=[fit])
        print(f"Wrote Emax bundle to {bundle}")
    return 0


def _ivive(
    input_json: Path,
    *,
    evidence_out: Path | None,
    study_dir: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.dmpk.ivive import (
        IviveRequest,
        predict_hepatic_clearance_well_stirred,
    )
    from pharmacoweave.integrations.pharmacology_evidence import (
        evidence_from_ivive,
        exchange_from_records,
        write_exchange,
    )

    if bundle is not None and study_dir is None:
        raise IngestError("--bundle requires --study so study provenance and validation are recorded.")
    if not input_json.is_file():
        raise IngestError(f"{input_json} does not exist.")
    try:
        raw = json.loads(input_json.read_text(encoding="utf-8"))
    except Exception as exc:
        raise IngestError(f"{input_json} is not valid JSON: {exc}") from exc

    if isinstance(raw, list):
        requests = [IviveRequest.model_validate(item) for item in raw]
    elif isinstance(raw, dict) and "requests" in raw:
        requests = [IviveRequest.model_validate(item) for item in raw["requests"]]
    else:
        requests = [IviveRequest.model_validate(raw)]

    predictions = []
    for req in requests:
        pred = predict_hepatic_clearance_well_stirred(
            req.clint_record,
            physiology=req.physiology,
            fu_plasma=req.fu_plasma,
            fu_inc=req.fu_inc,
            blood_to_plasma_ratio=req.blood_to_plasma_ratio,
            prediction_id=req.prediction_id,
        )
        predictions.append(pred)
        print(
            f"{pred.prediction_id}\tspecies={pred.species}\t"
            f"CLint_app={pred.clint_app_value:g} {pred.clint_app_unit}\t"
            f"CLH_plasma={pred.hepatic_plasma_clearance_l_h:.6g} L/h\t"
            f"CLH_plasma_kg={pred.hepatic_plasma_clearance_l_h_kg:.6g} L/h/kg\t"
            f"ER={pred.hepatic_extraction_ratio:.4g}"
        )

    selected_run = run_id or f"ivive-run-{predictions[0].compound_id}"
    if evidence_out is not None:
        evidence_records = [evidence_from_ivive(p, run_id=selected_run) for p in predictions]
        exchange = exchange_from_records(evidence_records)
        evidence_out.parent.mkdir(parents=True, exist_ok=True)
        write_exchange(exchange, evidence_out)
        print(f"Wrote PharmacologyEvidence to {evidence_out}")

    if bundle is not None and study_dir is not None:
        from pharmacoweave.core.study import PharmacologyStudy

        study = PharmacologyStudy.from_files(study_dir)
        study.export_bundle(bundle, run_id=selected_run, ivive_predictions=predictions)
        print(f"Wrote IVIVE bundle to {bundle}")
    return 0


def _audit_predictions(
    study_dir: Path,
    *,
    priors_path: Path | None,
    nca_config: Path | None,
    fit_model: str | None,
    estimates_path: Path | None,
    output: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.core.evidence import EvidenceRecord, EvidenceType
    from pharmacoweave.core.study import PharmacologyStudy
    from pharmacoweave.diagnostics.prediction_audit import audit_study_predictions
    from pharmacoweave.integrations.guards import load_exchange
    from pharmacoweave.integrations.pharmacology_evidence import collect_run_evidence

    study = PharmacologyStudy.from_files(study_dir)
    selected_run = run_id or f"{study.metadata.study_id}-audit"

    priors = (
        [r for r in load_exchange(priors_path).records if r.evidence_type is EvidenceType.PREDICTED]
        if priors_path is not None
        else [r for r in study.external_evidence if r.evidence_type is EvidenceType.PREDICTED]
    )

    estimates: list[EvidenceRecord] = []
    nca_results = None
    pk_fits: list[Any] | None = None
    if estimates_path is not None:
        est_exchange = load_exchange(estimates_path)
        estimates.extend(r for r in est_exchange.records if r.evidence_type is EvidenceType.ESTIMATED)
    if nca_config is not None:
        from pharmacoweave.pk_nca.runner import load_nca_config

        cfg = load_nca_config(nca_config)
        nca_results = study.run_nca(cfg)
        estimates.extend(collect_run_evidence(run_id=selected_run, nca_parameters=nca_results))
    if fit_model is not None:
        if fit_model == "one_compartment_iv_bolus":
            from pharmacoweave.pk_models.fit_iv import fit_iv_bolus

            pk_fits = list(fit_iv_bolus(study))
        elif fit_model == "one_compartment_iv_infusion":
            from pharmacoweave.pk_models.fit_infusion import fit_one_compartment_iv_infusion

            pk_fits = list(fit_one_compartment_iv_infusion(study))
        elif fit_model == "one_compartment_extravascular":
            from pharmacoweave.pk_models.fit_extravascular import fit_one_compartment_extravascular

            pk_fits = list(fit_one_compartment_extravascular(study))
        elif fit_model == "two_compartment_iv_bolus":
            from pharmacoweave.pk_models.fit_two_compartment import fit_two_compartment_iv_bolus

            pk_fits = list(fit_two_compartment_iv_bolus(study))
        elif fit_model == "three_compartment_iv_bolus":
            from pharmacoweave.pk_models.fit_three_compartment import fit_three_compartment_iv_bolus

            pk_fits = list(fit_three_compartment_iv_bolus(study))
        estimates.extend(collect_run_evidence(run_id=selected_run, pk_fits=pk_fits))

    estimates = [r for r in estimates if r.evidence_type is EvidenceType.ESTIMATED]

    audit_report = audit_study_predictions(study, estimates, priors=priors, audit_id=selected_run)

    for s in audit_report.summaries:
        print(
            f"{s.measurement_type}\tN={s.n_comparisons}\tAFE={s.afe:.4g}\tAAFE={s.aafe:.4g}\t"
            f"within_2_fold={s.fraction_within_2_fold * 100:.1f}%\t"
            f"within_3_fold={s.fraction_within_3_fold * 100:.1f}%\t"
            f"outliers={s.n_outlier_2_fold}"
        )
    for sk in audit_report.skipped:
        print(f"SKIPPED\t{sk.prior_evidence_id}\treason={sk.reason}")

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(audit_report.model_dump_json(indent=2) + "\n", encoding="utf-8")
        print(f"Wrote PredictionAuditReport to {output}")

    if bundle is not None:
        study.export_bundle(
            bundle,
            run_id=selected_run,
            nca_parameters=nca_results,
            fits=pk_fits,
            prediction_audit=audit_report,
        )
        print(f"Wrote audit bundle to {bundle}")
    return 0


def _verify_bundle(bundle_dir: Path, *, as_json: bool) -> int:
    from pharmacoweave.core.bundle_verify import verify_bundle

    result = verify_bundle(bundle_dir)
    if as_json:
        print(result.model_dump_json(indent=2))
    else:
        for line in result.summary_lines():
            print(line)
    return 0 if result.is_valid else 1


def _diff(bundle_a: Path, bundle_b: Path, *, as_json: bool, rel_tol: float) -> int:
    from pharmacoweave.core.bundle_diff import diff_bundles

    result = diff_bundles(bundle_a, bundle_b, rel_tol=rel_tol)
    if as_json:
        print(result.model_dump_json(indent=2))
    else:
        for line in result.summary_lines():
            print(line)
    return 0 if result.identical else 1


def _ingest(
    mapping_path: Path,
    *,
    output_dir: Path,
    file_format: str,
    strict: bool,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.ingestion.mapping import ingest_study_from_mapping
    from pharmacoweave.ingestion.tabular import write_study_directory

    mode = "strict" if strict else "research"
    study = ingest_study_from_mapping(mapping_path, mode=mode, validate_schema=True)
    write_study_directory(study, output_dir, file_format=file_format)
    print(f"Wrote canonical study directory to {output_dir}")

    if bundle is not None:
        selected_run = run_id or f"{study.metadata.study_id}-ingest"
        study.export_bundle(bundle, run_id=selected_run, mode=mode)
        print(f"Wrote result bundle to {bundle}")
    return 0


def _load_simulation_request(input_json: Path):
    import json

    from pydantic import ValidationError

    from pharmacoweave.simulation.exposure import ExposureSimulationRequest

    try:
        raw = json.loads(input_json.read_text(encoding="utf-8"))
    except Exception as exc:
        raise IngestError(f"{input_json} is not valid JSON: {exc}") from exc
    try:
        return ExposureSimulationRequest.model_validate(raw)
    except ValidationError as exc:
        raise IngestError(f"{input_json} is not an exposure simulation request: {exc}") from exc


def _accepted_record(
    request_json: Path,
    *,
    output: Path | None,
    study_dir: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    import json

    from pharmacoweave.accepted_records import evaluate_accepted_record
    from pharmacoweave.core.bundle import _json_ready

    try:
        payload = json.loads(request_json.read_text(encoding="utf-8"))
    except Exception as exc:
        raise IngestError(f"{request_json} is not valid JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise IngestError(f"{request_json} must contain a JSON object.")
    result = evaluate_accepted_record(payload)
    ready = _json_ready(result)
    print(f"kind={payload.get('kind')}\tdosing_recommendation=false")
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(ready, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"Wrote accepted record to {output}")
    if bundle is not None and study_dir is not None:
        from pharmacoweave.core.study import PharmacologyStudy

        study = PharmacologyStudy.from_files(study_dir)
        selected_run = run_id or "accepted-record"
        study.export_bundle(bundle, run_id=selected_run, accepted_record=result)
        print(f"Wrote accepted record bundle to {bundle}")
    elif bundle is not None:
        raise IngestError("accepted-record --bundle requires --study.")
    return 0


def _simulate(
    input_json: Path,
    *,
    output: Path | None,
    study_dir: Path | None,
    bundle: Path | None,
    run_id: str | None,
) -> int:
    from pharmacoweave.simulation.exposure import simulate_exposure

    request = _load_simulation_request(input_json)
    result = simulate_exposure(request)
    exposure = (
        f"AUC_inf={result.auc_inf:.6g} {result.auc_unit}"
        if result.auc_inf is not None
        else f"AUC_tau={result.auc_tau:.6g} {result.auc_unit}"
    )
    print(
        f"{result.simulation_id}\tmodel={result.model_id}\tregimen={result.regimen}\t"
        f"Cmax={result.cmax:.6g} {result.concentration_unit}\t"
        f"{exposure}\t"
        f"evidence={result.evidence_type}\tdosing_recommendation=false"
    )
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(result.model_dump_json(indent=2), encoding="utf-8")
        print(f"Wrote exposure simulation to {output}")
    if bundle is not None and study_dir is not None:
        from pharmacoweave.core.study import PharmacologyStudy

        study = PharmacologyStudy.from_files(study_dir)
        selected_run = run_id or f"{request.simulation_id}-simulation"
        study.export_bundle(bundle, run_id=selected_run, exposure_simulation=result)
        print(f"Wrote simulation bundle to {bundle}")
    elif bundle is not None:
        raise IngestError("simulate --bundle requires --study.")
    return 0


def _engine_run(input_json: Path, *, engine: str, output: Path | None) -> int:
    from pharmacoweave.integrations.external_backends import run_backend

    request = _load_simulation_request(input_json)
    response = run_backend(engine, request)
    result = response.result
    exposure = (
        f"AUC_inf={result.auc_inf:.6g} {result.auc_unit}"
        if result.auc_inf is not None
        else f"AUC_tau={result.auc_tau:.6g} {result.auc_unit}"
    )
    print(
        f"{response.engine_name}\t{response.license_spdx}\tcopyleft={str(response.copyleft).lower()}\t"
        f"Cmax={result.cmax:.6g} {result.concentration_unit}\t"
        f"{exposure}\tdosing_recommendation=false"
    )
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(response.model_dump_json(indent=2), encoding="utf-8")
        print(f"Wrote engine response to {output}")
    return 0


def _engine_slot(
    request_file: Path,
    *,
    executable: Path,
    engine_version: str,
    license_file: Path,
    output: Path | None,
) -> int:
    from pharmacoweave.integrations.external_backends import run_caller_engine

    if not request_file.is_file():
        raise IngestError(f"Caller engine request {request_file} does not exist.")
    if not license_file.is_file():
        raise IngestError(f"Caller engine license file {license_file} does not exist.")
    recorded = run_caller_engine(
        executable,
        request_file.read_text(encoding="utf-8"),
        engine_version=engine_version,
        license_text=license_file.read_text(encoding="utf-8"),
    )
    print(
        f"executable_sha256={recorded.executable_sha256}\tengine_version={recorded.engine_version}\t"
        f"response_sha256={recorded.response_sha256}\tdosing_recommendation=false"
    )
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(recorded.model_dump_json(indent=2), encoding="utf-8")
        print(f"Wrote caller engine run to {output}")
    return 0


def _capabilities() -> int:
    for item in MODULE_STATUS:
        print(f"{item.name:16} {item.phase:4} {item.execution:10} {item.summary}")
    return 0


def _export_schemas(directory: Path) -> int:
    from pharmacoweave.schema_export import export_schemas

    export_schemas(directory)
    print(f"Wrote JSON Schemas to {directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
