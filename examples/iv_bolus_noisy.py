"""Noisy single-dose IV bolus used to recover CL and V.

One subject, dose 1 mg, CL = 0.40 L/h, V = 2.0 L. Quantified concentrations
were drawn once from C(t) * exp(0.10 * z) with ``random.Random(20261030)``
and stored to 8 decimal places. The 24 h sample is below 0.01 µg/mL and is
BLQ with no numeric value. Body weight is recorded; CL and V are absolute,
not per kilogram.

The relative recovery tolerance for this draw is 10% on both CL and V.
"""

from datetime import UTC, datetime, timedelta
from pathlib import Path

from pharmacoweave.core.entities import (
    Analyte,
    AnalyteForm,
    AnalyteRole,
    Assay,
    Biosample,
    Compound,
    DataClass,
    DoseEvent,
    Observation,
    ObservationStatus,
    ObservationType,
    QcStatus,
    Quantity,
    Route,
    Sex,
    StudyMetadata,
    StudyType,
    Subject,
    Timebase,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.ingestion.tabular import write_study_directory

DOSE_TIME = datetime(2024, 3, 2, 9, 0, tzinfo=UTC)
TRUE_CLEARANCE_L_H = 0.40
TRUE_VOLUME_L = 2.0
DOSE_MG = 1.0
BODY_KG = 0.25
LLOQ_UG_PER_ML = 0.01
RECOVERY_RELATIVE_TOLERANCE = 0.10
# time in minutes, concentration in ug/mL (numerically equal to mg/L).
QUANTIFIED = (
    (15, 0.42056065),
    (30, 0.43641680),
    (60, 0.41693394),
    (120, 0.32831068),
    (240, 0.23287745),
    (360, 0.13239253),
    (480, 0.10034908),
    (720, 0.04909356),
)
EXAMPLE_DIR = Path(__file__).resolve().parent / "data" / "iv_bolus_noisy"


def build_study() -> PharmacologyStudy:
    """Build the noisy synthetic study in memory."""

    metadata = StudyMetadata(
        study_id="PW-SYNTH-IV-NOISY-001",
        protocol_id="SYNTH-PROTOCOL-NOISY-001",
        title="Synthetic noisy single-dose IV bolus plasma PK",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Rattus norvegicus",
        ncbi_taxon_id=10116,
        design_version="1",
        description=(
            "Synthetic plasma concentrations for a one-compartment IV bolus with "
            "multiplicative noise. Dose = 1 mg, CL = 0.40 L/h, V = 2.0 L. "
            "Quantified values were drawn once from C(t) * exp(0.10 * z) using "
            "random.Random(20261030) and rounded to 8 decimal places. The 24 h "
            "sample is below the 0.01 ug/mL LLOQ and is stored as BLQ. This is "
            "not an animal experiment and not a dosing recommendation."
        ),
        timebases=[
            Timebase(timebase_id="STUDY-CLOCK", description="Synthetic protocol clock, UTC."),
        ],
    )
    subject = Subject(
        subject_id="SUBJ-01",
        study_id=metadata.study_id,
        group_id="IV-NOISY",
        species=metadata.species,
        sex=Sex.FEMALE,
        body_weight=Quantity(value=BODY_KG, unit="kg"),
    )
    compound = Compound(compound_id="PW-0001", name="Example compound PW-0001")
    analyte = Analyte(
        analyte_id="PW-0001-PARENT",
        compound_id=compound.compound_id,
        name="PW-0001 parent",
        analyte_role=AnalyteRole.PARENT,
        matrix="plasma",
    )
    assay = Assay(
        assay_id="ASSAY-LCMS-01",
        method="LC-MS/MS",
        analyte_id=analyte.analyte_id,
        matrix="plasma",
        calibration_id="CAL-NOISY-001",
        lloq=LLOQ_UG_PER_ML,
        uloq=10.0,
        limit_unit="ug/mL",
        qc_status=QcStatus.ACCEPTED,
    )
    dose = DoseEvent(
        dose_event_id="DOSE-SUBJ-01",
        subject_id=subject.subject_id,
        compound_id=compound.compound_id,
        absolute_time=DOSE_TIME,
        dose_value=DOSE_MG,
        dose_unit="mg",
        route=Route.IV_BOLUS,
        formulation="saline solution",
        body_weight_normalized_value=DOSE_MG / BODY_KG,
        body_weight_normalized_unit="mg/kg",
    )
    samples: list[Biosample] = []
    observations: list[Observation] = []
    for minute in (-5, *(time for time, _value in QUANTIFIED), 1440):
        sample_id = "SUBJ-01-PRE" if minute < 0 else f"SUBJ-01-T{minute:04d}"
        samples.append(
            Biosample(
                sample_id=sample_id,
                subject_id=subject.subject_id,
                collection_time_absolute=DOSE_TIME + timedelta(minutes=minute),
                time_since_dose=float(minute),
                time_since_dose_unit="min",
                dose_event_id=dose.dose_event_id,
                matrix="plasma",
                collection_site="venous",
                timebase_id="STUDY-CLOCK",
                is_predose=minute < 0,
            )
        )
        observations.append(_observation(sample_id, minute, assay.assay_id, analyte.analyte_id))
    return PharmacologyStudy(
        metadata=metadata,
        subjects=[subject],
        compounds=[compound],
        analytes=[analyte],
        assays=[assay],
        dose_events=[dose],
        biosamples=samples,
        observations=observations,
    )


def _observation(sample_id: str, minute: int, assay_id: str, analyte_id: str) -> Observation:
    common = {
        "observation_id": f"OBS-{sample_id}",
        "sample_id": sample_id,
        "analyte_id": analyte_id,
        "observation_type": ObservationType.CONCENTRATION,
        "assay_id": assay_id,
        "analyte_form": AnalyteForm.PARENT,
        "assay_batch_id": "BATCH-NOISY-01",
        "unit": "ug/mL",
    }
    quantified = dict(QUANTIFIED)
    if minute < 0 or minute not in quantified:
        return Observation(
            **common,
            status=ObservationStatus.BLQ,
            raw_value="<0.01",
            reported_value="<LLOQ",
        )
    text = f"{quantified[minute]:.8f}"
    return Observation(
        **common,
        status=ObservationStatus.QUANTIFIED,
        raw_value=text,
        reported_value=text,
        numeric_value=quantified[minute],
    )


def main() -> None:
    """Write the example study directory."""

    write_study_directory(build_study(), EXAMPLE_DIR)
    print(f"Wrote {EXAMPLE_DIR}")


if __name__ == "__main__":
    main()
