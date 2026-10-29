"""Noise-free single-dose IV bolus study used as the Phase 0 fixture.

Two subjects share CL = 0.20 L/h/kg and V = 1.0 L/kg. Concentrations follow
the one-compartment IV bolus equation. The 24 h sample is below the LLOQ and
is stored as BLQ, with no numeric substitute.
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
from pharmacoweave.core.evidence import IdentityLink, ProjectName, ReviewStatus
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.ingestion.tabular import write_study_directory
from pharmacoweave.pk_models.analytical import iv_bolus_concentration

DOSE_TIME = datetime(2024, 1, 15, 9, 0, tzinfo=UTC)
CLEARANCE_L_PER_H_PER_KG = 0.20
VOLUME_L_PER_KG = 1.0
DOSE_MG_PER_KG = 1.0
LLOQ_UG_PER_ML = 0.01
ULOQ_UG_PER_ML = 10.0
POSTDOSE_MINUTES = (5, 15, 30, 60, 120, 240, 480, 720, 1440)
SUBJECTS = (
    ("SUBJ-01", Sex.FEMALE, 0.25),
    ("SUBJ-02", Sex.MALE, 0.30),
)
EXAMPLE_DIR = Path(__file__).resolve().parent / "data" / "iv_bolus_single_dose"


def build_study() -> PharmacologyStudy:
    """Build the synthetic study in memory."""

    metadata = StudyMetadata(
        study_id="PW-SYNTH-IV-001",
        protocol_id="SYNTH-PROTOCOL-001",
        title="Synthetic single-dose IV bolus plasma PK",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Rattus norvegicus",
        ncbi_taxon_id=10116,
        design_version="1",
        description=(
            "Noise-free synthetic plasma concentrations for a one-compartment IV bolus. "
            "Both subjects use CL = 0.20 L/h/kg and V = 1.0 L/kg. The 24 h sample is below "
            "the 0.01 ug/mL LLOQ and is stored as BLQ. This is not an animal experiment "
            "and not a dosing recommendation."
        ),
        timebases=[
            Timebase(timebase_id="STUDY-CLOCK", description="Synthetic protocol clock, UTC."),
        ],
    )
    subjects = [
        Subject(
            subject_id=subject_id,
            study_id=metadata.study_id,
            group_id="IV-1mgkg",
            species=metadata.species,
            sex=sex,
            body_weight=Quantity(value=body_kg, unit="kg"),
        )
        for subject_id, sex, body_kg in SUBJECTS
    ]
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
        calibration_id="CAL-001",
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
                route=Route.IV_BOLUS,
                formulation="saline solution",
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
                _observation(subject_id, sample_id, minute, body_kg, assay.assay_id, analyte.analyte_id)
            )
    link = IdentityLink(
        link_id="LINK-EW-SUBJ-01",
        source_project=ProjectName.EPHYSWEAVE,
        source_study_id="EW-SYNTH-001",
        source_subject_id="EW-RAT-A",
        canonical_subject_id="SUBJ-01",
        match_method="synthetic_protocol",
        review_status=ReviewStatus.ACCEPTED,
        notes=(
            "Synthetic cross-reference. Accepted because this protocol defined the map. "
            "Matching text was not treated as identity."
        ),
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
        identity_links=[link],
    )


def expected_concentration_ug_ml(body_kg: float, minute: int) -> float:
    """Return the noise-free concentration in ug/mL."""

    return iv_bolus_concentration(
        dose=body_kg * DOSE_MG_PER_KG,
        dose_unit="mg",
        volume=body_kg * VOLUME_L_PER_KG,
        volume_unit="L",
        clearance=body_kg * CLEARANCE_L_PER_H_PER_KG,
        clearance_unit="L/h",
        time=float(minute),
        time_unit="min",
        concentration_unit="ug/mL",
    )


def _observation(
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
        "assay_batch_id": "BATCH-01",
        "unit": "ug/mL",
    }
    if minute < 0:
        return Observation(
            **common,
            status=ObservationStatus.BLQ,
            raw_value="<0.01",
            reported_value="<LLOQ",
        )
    concentration = expected_concentration_ug_ml(body_kg, minute)
    if concentration < LLOQ_UG_PER_ML:
        return Observation(
            **common,
            status=ObservationStatus.BLQ,
            raw_value="<0.01",
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
    """Write the example study directory."""

    write_study_directory(build_study(), EXAMPLE_DIR)
    print(f"Wrote {EXAMPLE_DIR}")


if __name__ == "__main__":
    main()
