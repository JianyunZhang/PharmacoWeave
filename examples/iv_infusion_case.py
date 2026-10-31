"""Synthetic single-dose constant-rate IV infusion study.

Two subjects receive 2.0 mg/kg over a 60 min constant-rate IV infusion with
CL = 0.25 L/h/kg and V = 1.0 L/kg. Predose (-5 min) and 24 h (1440 min) samples
are below the 0.02 ug/mL LLOQ and are stored as BLQ without numeric substitution.
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
from pharmacoweave.pk_models.analytical import iv_infusion_concentration

DOSE_TIME = datetime(2024, 3, 5, 9, 0, tzinfo=UTC)
INFUSION_DURATION_MIN = 60.0
CLEARANCE_L_PER_H_PER_KG = 0.25
VOLUME_L_PER_KG = 1.0
DOSE_MG_PER_KG = 2.0
MOLECULAR_WEIGHT_G_PER_MOL = 250.0
LLOQ_UG_PER_ML = 0.02
ULOQ_UG_PER_ML = 10.0
POSTDOSE_MINUTES = (15, 30, 60, 90, 120, 240, 480, 720, 1440)
SUBJECTS = (
    ("SUBJ-01", Sex.FEMALE, 0.25),
    ("SUBJ-02", Sex.MALE, 0.30),
)
INFUSION_EXAMPLE_DIR = Path(__file__).resolve().parent / "data" / "iv_infusion_single_dose"


def build_iv_infusion_study() -> PharmacologyStudy:
    """Build the noise-free single-dose constant-rate IV infusion study in memory."""

    metadata = StudyMetadata(
        study_id="PW-SYNTH-INF-001",
        protocol_id="SYNTH-PROTOCOL-INF-001",
        title="Synthetic single-dose constant-rate IV infusion plasma PK",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Rattus norvegicus",
        ncbi_taxon_id=10116,
        design_version="1",
        description=(
            "Noise-free synthetic plasma concentrations for a 60 min constant-rate IV infusion "
            "(CL = 0.25 L/h/kg, V = 1.0 L/kg, dose = 2.0 mg/kg). Predose and 24 h samples are "
            "below the 0.02 ug/mL LLOQ and stored as BLQ. Research-use synthetic fixture only."
        ),
        timebases=[Timebase(timebase_id="STUDY-CLOCK", description="Synthetic protocol clock, UTC.")],
    )
    subjects = [
        Subject(
            subject_id=subject_id,
            study_id=metadata.study_id,
            group_id="IVINF-2mgkg",
            species=metadata.species,
            sex=sex,
            body_weight=Quantity(value=body_kg, unit="kg"),
        )
        for subject_id, sex, body_kg in SUBJECTS
    ]
    compound = Compound(
        compound_id="PW-0001",
        name="Example compound PW-0001",
        molecular_weight=MOLECULAR_WEIGHT_G_PER_MOL,
        molecular_weight_unit="g/mol",
    )
    analyte = Analyte(
        analyte_id="PW-0001-PARENT",
        compound_id=compound.compound_id,
        name="PW-0001 parent",
        analyte_role=AnalyteRole.PARENT,
        matrix="plasma",
    )
    assay = Assay(
        assay_id="ASSAY-LCMS-INF-01",
        method="LC-MS/MS",
        analyte_id=analyte.analyte_id,
        matrix="plasma",
        calibration_id="CAL-INF-001",
        lloq=LLOQ_UG_PER_ML,
        uloq=ULOQ_UG_PER_ML,
        limit_unit="ug/mL",
        qc_status=QcStatus.ACCEPTED,
    )
    doses: list[DoseEvent] = []
    samples: list[Biosample] = []
    observations: list[Observation] = []
    for subject_id, _sex, body_kg in SUBJECTS:
        doses.append(
            DoseEvent(
                dose_event_id=f"DOSE-{subject_id}",
                subject_id=subject_id,
                compound_id=compound.compound_id,
                absolute_time=DOSE_TIME,
                dose_value=body_kg * DOSE_MG_PER_KG,
                dose_unit="mg",
                route=Route.IV_INFUSION,
                formulation="saline infusion",
                infusion_duration=INFUSION_DURATION_MIN,
                infusion_duration_unit="min",
                body_weight_normalized_value=DOSE_MG_PER_KG,
                body_weight_normalized_unit="mg/kg",
            )
        )
        for minute in (-5, *POSTDOSE_MINUTES):
            sample_id = f"{subject_id}-PRE" if minute < 0 else f"{subject_id}-T{minute:04d}"
            samples.append(
                Biosample(
                    sample_id=sample_id,
                    subject_id=subject_id,
                    collection_time_absolute=DOSE_TIME + timedelta(minutes=minute),
                    time_since_dose=float(minute),
                    time_since_dose_unit="min",
                    dose_event_id=f"DOSE-{subject_id}",
                    matrix="plasma",
                    collection_site="venous",
                    timebase_id="STUDY-CLOCK",
                    is_predose=minute < 0,
                )
            )
            observations.append(
                _infusion_observation(subject_id, sample_id, minute, body_kg, assay.assay_id, analyte.analyte_id)
            )
    return PharmacologyStudy(
        metadata=metadata,
        subjects=subjects,
        compounds=[compound],
        analytes=[analyte],
        assays=[assay],
        dose_events=doses,
        biosamples=samples,
        observations=observations,
    )


def expected_infusion_concentration_ug_ml(body_kg: float, minute: int) -> float:
    """Return the noise-free IV infusion concentration in ug/mL."""

    return iv_infusion_concentration(
        dose=body_kg * DOSE_MG_PER_KG,
        dose_unit="mg",
        infusion_duration=INFUSION_DURATION_MIN,
        infusion_duration_unit="min",
        volume=body_kg * VOLUME_L_PER_KG,
        volume_unit="L",
        clearance=body_kg * CLEARANCE_L_PER_H_PER_KG,
        clearance_unit="L/h",
        time=float(minute),
        time_unit="min",
        concentration_unit="ug/mL",
    )


def _infusion_observation(
    subject_id: str,
    sample_id: str,
    minute: int,
    body_kg: float,
    assay_id: str,
    analyte_id: str,
) -> Observation:
    common = {
        "observation_id": f"OBS-{sample_id}",
        "sample_id": sample_id,
        "analyte_id": analyte_id,
        "observation_type": ObservationType.CONCENTRATION,
        "assay_id": assay_id,
        "analyte_form": AnalyteForm.PARENT,
        "assay_batch_id": "BATCH-INF-01",
        "unit": "ug/mL",
    }
    if minute < 0:
        return Observation(
            **common,
            status=ObservationStatus.BLQ,
            raw_value="<0.02",
            reported_value="<LLOQ",
        )
    concentration = expected_infusion_concentration_ug_ml(body_kg, minute)
    if concentration < LLOQ_UG_PER_ML:
        return Observation(
            **common,
            status=ObservationStatus.BLQ,
            raw_value="<0.02",
            reported_value="<LLOQ",
        )
    rounded = round(concentration, 6)
    text = f"{rounded:.6f}"
    return Observation(
        **common,
        status=ObservationStatus.QUANTIFIED,
        raw_value=text,
        reported_value=text,
        numeric_value=rounded,
    )


def main() -> None:
    """Write the synthetic single-dose IV infusion study directory."""

    write_study_directory(build_iv_infusion_study(), INFUSION_EXAMPLE_DIR)
    print(f"Wrote {INFUSION_EXAMPLE_DIR}")


if __name__ == "__main__":
    main()
