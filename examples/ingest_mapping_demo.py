"""Demonstration of declarative ingest mapping from heterogeneous lab tables.

Reads messy lab CSV files (animals.csv, dosing.csv, pk_concentrations.csv)
according to ingest_mapping.yaml, preserves raw SHA-256 file hashes,
validates the study relational model, and writes a ResultBundle.
"""

from pathlib import Path

from pharmacoweave.core.bundle_verify import verify_bundle
from pharmacoweave.core.study import PharmacologyStudy


def run_demo() -> None:
    demo_dir = Path(__file__).resolve().parent / "data" / "ingest_mapping_demo"
    mapping_file = demo_dir / "ingest_mapping.yaml"

    print("=== PharmacoWeave 0.1.0: Declarative Ingest Mapping Demo ===")
    print(f"Reading mapping configuration from: {mapping_file}")

    # Load study via classmethod
    study = PharmacologyStudy.from_mapping(mapping_file)
    print(f"\nIngested Study: {study.metadata.study_id} ({study.metadata.title})")
    print(f"  Subjects:        {len(study.subjects)} ([{', '.join(s.subject_id for s in study.subjects)}])")
    print(f"  Dosing Events:   {len(study.dose_events)} ([{', '.join(d.dose_event_id for d in study.dose_events)}])")
    print(f"  Biosamples:      {len(study.biosamples)}")
    n_blq = sum(1 for o in study.observations if o.status == "BLQ")
    print(f"  Observations:    {len(study.observations)} ({n_blq} BLQ)")

    print("\nCryptographic Raw Source File Provenance:")
    for sf in study.source_files:
        print(f"  [{sf.role:16}] {sf.relative_path:24} SHA-256: {sf.sha256}")

    # Validate study
    report = study.check()
    print(f"\nValidation Report: ok={report.ok}, issues={len(report.issues)}")

    # Export a bundle to a temporary location and verify it
    out_bundle = demo_dir / "output_bundle"
    if out_bundle.exists():
        import shutil
        shutil.rmtree(out_bundle)

    manifest = study.export_bundle(out_bundle, run_id="ingest-demo-run-01")
    print(f"\nExported ResultBundle to {out_bundle} (files: {len(manifest.files)})")

    # Offline verify the written bundle
    verification = verify_bundle(out_bundle)
    print(f"Bundle Verification Result: valid={verification.is_valid}, files_checked={verification.files_checked}")

    # Clean up output bundle
    import shutil
    shutil.rmtree(out_bundle)
    print("\nDemonstration complete.")


if __name__ == "__main__":
    run_demo()
