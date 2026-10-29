"""Command line for validation, the locked IV bolus analysis, and capabilities."""

import argparse
from pathlib import Path

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

    analyze = commands.add_parser("analyze", help="Run locked IV bolus NCA and the one-compartment fit")
    analyze.add_argument("study_dir", type=Path)
    analyze.add_argument("--config", type=Path, required=True)
    analyze.add_argument("--bundle", type=Path, required=True)
    analyze.add_argument("--run-id")
    analyze.add_argument("--strict", action="store_true")

    commands.add_parser("capabilities", help="Show which modules can execute")
    export = commands.add_parser("export-schemas", help="Write JSON Schema files from the data contracts")
    export.add_argument("--directory", type=Path, default=Path("schemas"))

    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _validate(args.study_dir, strict=args.strict, bundle=args.bundle, run_id=args.run_id)
        if args.command == "analyze":
            return _analyze(
                args.study_dir,
                config=args.config,
                bundle=args.bundle,
                run_id=args.run_id,
                strict=args.strict,
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


def _analyze(study_dir: Path, *, config: Path, bundle: Path, run_id: str | None, strict: bool) -> int:
    from pharmacoweave.core.study import PharmacologyStudy

    mode = "strict" if strict else "research"
    study = PharmacologyStudy.from_files(study_dir, mode=mode)
    selected_run = run_id or f"{study.metadata.study_id}-analysis"
    manifest = study.export_analysis_bundle(bundle, run_id=selected_run, nca_config=config, mode=mode)
    print(f"analysis_executed={manifest.analysis_executed}")
    print(f"Wrote analysis bundle to {bundle}")
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
