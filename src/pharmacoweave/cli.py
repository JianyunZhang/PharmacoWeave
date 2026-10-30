"""Command line for validation, NCA, PK/PD fitting, schema export, and capabilities."""

import argparse
import json
from pathlib import Path
from typing import Literal

from pharmacoweave.core.roadmap import MODULE_STATUS
from pharmacoweave.errors import IngestError, PharmacoWeaveError, ValidationFailed
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

    nca = commands.add_parser("nca", help="Run single-dose IV bolus NCA")
    nca.add_argument("study_dir", type=Path)
    nca.add_argument("--config", type=Path, required=True)
    nca.add_argument(
        "--with-fit",
        action="store_true",
        help="Also run log-linear one-compartment IV bolus fit alongside NCA when exporting a bundle",
    )
    nca.add_argument(
        "--bundle",
        type=Path,
        help="Write a bundle with NCA tables and concentration-time figures",
    )
    nca.add_argument("--run-id")

    fit = commands.add_parser("fit", help="Fit a one-compartment IV bolus for CL and V")
    fit.add_argument("study_dir", type=Path)
    fit.add_argument("--method", choices=("log_linear_ols", "weighted_nls"), default="log_linear_ols")
    fit.add_argument("--weight", choices=("1/y", "1/y^2"), default="1/y^2")
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

    commands.add_parser("capabilities", help="Show which modules can execute")
    export = commands.add_parser("export-schemas", help="Write JSON Schema files from the data contracts")
    export.add_argument("--directory", type=Path, default=Path("schemas"))

    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _validate(args.study_dir, strict=args.strict, bundle=args.bundle, run_id=args.run_id)
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
                method=args.method,
                weight=args.weight,
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
    from pharmacoweave.core.study import PharmacologyStudy

    study = PharmacologyStudy.from_files(study_dir)
    results = study.run_nca(config)
    for row in results:
        print(
            f"{row.subject_id}\tAUC_inf={row.auc_inf}\tCL={row.clearance}\thalf_life={row.half_life}\tVz={row.volume}"
        )
    if bundle is not None:
        fits = None
        if with_fit:
            from pharmacoweave.pk_models.fit_iv import fit_iv_bolus

            fits = fit_iv_bolus(study)
        selected_run = run_id or f"{study.metadata.study_id}-nca"
        study.export_bundle(bundle, run_id=selected_run, nca_parameters=results, fits=fits)
        print(f"Wrote NCA bundle to {bundle}")
    return 0


def _fit(study_dir: Path, *, method: str, weight: str, bundle: Path | None, run_id: str | None) -> int:
    from pharmacoweave.core.study import PharmacologyStudy
    from pharmacoweave.pk_models.fit_iv import fit_iv_bolus
    from pharmacoweave.pk_models.wnls import fit_iv_bolus_wnls

    study = PharmacologyStudy.from_files(study_dir)
    fits = fit_iv_bolus(study) if method == "log_linear_ols" else fit_iv_bolus_wnls(study, weight=weight)
    for item in fits:
        values = {estimate.parameter: estimate.value for estimate in item.estimates}
        clearance = "" if "CL" not in values else f"{values['CL']:.6g}"
        volume = "" if "V" not in values else f"{values['V']:.6g}"
        print(
            f"{item.subject_id}\tmethod={item.fit.estimation_method}\t"
            f"status={item.diagnostic.status.value}\tCL={clearance}\tV={volume}"
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
