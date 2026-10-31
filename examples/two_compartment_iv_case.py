"""Synthetic two-compartment IV bolus benchmark cases (quiet and noisy).

True central/peripheral parameters for SUBJ-01 (Dose = 10.0 mg):
CL = 0.50 L/h, V1 = 1.0 L, Q = 1.50 L/h, V2 = 3.0 L.
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
from pharmacoweave.core.profiles import ConcentrationProfile, QuantifiedPoint
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.ingestion.tabular import write_study_directory
from pharmacoweave.pk_models.analytical import two_compartment_iv_bolus_concentration

DOSE_TIME = datetime(2024, 4, 1, 9, 0, tzinfo=UTC)
DOSE_MG = 10.0
TRUE_CL_L_H = 0.50
TRUE_V1_L = 1.00
TRUE_Q_L_H = 1.50
TRUE_V2_L = 3.00
LLOQ_UG_PER_ML = 0.05
ULOQ_UG_PER_ML = 25.0
POSTDOSE_MINUTES = (5, 10, 20, 40, 60, 120, 240, 360, 480, 720, 960, 1440, 2880)
TWO_COMP_EXAMPLE_DIR = Path(__file__).resolve().parent / "data" / "two_compartment_iv_bolus"

_NOISY_MULTIPLIERS = (1.03, 0.97, 1.04, 0.96, 1.02, 0.98, 1.03, 0.97, 1.02, 0.96, 1.03, 0.98)


def build_two_compartment_study() -> PharmacologyStudy:
    """Build the noise-free two-compartment IV bolus study in memory."""

    metadata = StudyMetadata(
        study_id="PW-SYNTH-2C-IV-001",
        protocol_id="SYNTH-PROTOCOL-2C-001",
        title="Synthetic two-compartment IV bolus plasma PK",
        study_type=StudyType.PK,
        data_class=DataClass.SYNTHETIC,
        species="Canis lupus familiaris",
        ncbi_taxon_id=9615,
        design_version="1",
        description=(
            "Noise-free synthetic biexponential IV bolus plasma concentrations "
            "(Dose = 10 mg, CL = 0.50 L/h, V1 = 1.0 L, Q = 1.50 L/h, V2 = 3.0 L). "
            "The 48 h sample is below the 0.05 ug/mL LLOQ and stored as BLQ. "
            "Research-use synthetic fixture only."
        ),
        timebases=[Timebase(timebase_id="STUDY-CLOCK", description="Synthetic protocol clock, UTC.")],
    )
    subject = Subject(
        subject_id="SUBJ-01",
        study_id=metadata.study_id,
        group_id="IV-10mg",
        species=metadata.species,
        sex=Sex.FEMALE,
        body_weight=Quantity(value=10.0, unit="kg"),
    )
    compound = Compound(
        compound_id="PW-0002",
        name="Example biexponential compound PW-0002",
        molecular_weight=400.0,
        molecular_weight_unit="g/mol",
    )
    analyte = Analyte(
        analyte_id="PW-0002-PARENT",
        compound_id=compound.compound_id,
        name="PW-0002 parent",
        analyte_role=AnalyteRole.PARENT,
        matrix="plasma",
    )
    assay = Assay(
        assay_id="ASSAY-LCMS-2C-01",
        method="LC-MS/MS",
        analyte_id=analyte.analyte_id,
        matrix="plasma",
        calibration_id="CAL-2C-001",
        lloq=LLOQ_UG_PER_ML,
        uloq=ULOQ_UG_PER_ML,
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
        formulation="saline bolus",
        body_weight_normalized_value=1.0,
        body_weight_normalized_unit="mg/kg",
    )
    samples: list[Biosample] = []
    observations: list[Observation] = []
    for minute in (-5, *POSTDOSE_MINUTES):
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
        common = {
            "observation_id": f"OBS-{sample_id}",
            "sample_id": sample_id,
            "analyte_id": analyte.analyte_id,
            "observation_type": ObservationType.CONCENTRATION,
            "assay_id": assay.assay_id,
            "analyte_form": AnalyteForm.PARENT,
            "assay_batch_id": "BATCH-2C-01",
            "unit": "ug/mL",
        }
        if minute < 0:
            observations.append(
                Observation(**common, status=ObservationStatus.BLQ, raw_value="<0.05", reported_value="<LLOQ")
            )
            continue
        conc = two_compartment_iv_bolus_concentration(
            dose=DOSE_MG,
            dose_unit="mg",
            clearance=TRUE_CL_L_H,
            clearance_unit="L/h",
            central_volume=TRUE_V1_L,
            central_volume_unit="L",
            intercompartmental_clearance=TRUE_Q_L_H,
            intercompartmental_clearance_unit="L/h",
            peripheral_volume=TRUE_V2_L,
            peripheral_volume_unit="L",
            time=float(minute),
            time_unit="min",
            concentration_unit="ug/mL",
        )
        if conc < LLOQ_UG_PER_ML:
            observations.append(
                Observation(**common, status=ObservationStatus.BLQ, raw_value="<0.05", reported_value="<LLOQ")
            )
        else:
            rounded = round(conc, 6)
            text = f"{rounded:.6f}"
            observations.append(
                Observation(
                    **common,
                    status=ObservationStatus.QUANTIFIED,
                    raw_value=text,
                    reported_value=text,
                    numeric_value=rounded,
                )
            )
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


def build_noisy_two_compartment_profile() -> ConcentrationProfile:
    """Build a deterministic noisy biexponential IV bolus profile."""

    quantified_minutes = POSTDOSE_MINUTES[:-1]
    points = tuple(
        QuantifiedPoint(
            observation_id=f"OBS-2C-NOISY-{idx:02d}",
            time_h=minute / 60.0,
            concentration_mg_l=round(
                mult
                * two_compartment_iv_bolus_concentration(
                    dose=DOSE_MG,
                    dose_unit="mg",
                    clearance=TRUE_CL_L_H,
                    clearance_unit="L/h",
                    central_volume=TRUE_V1_L,
                    central_volume_unit="L",
                    intercompartmental_clearance=TRUE_Q_L_H,
                    intercompartmental_clearance_unit="L/h",
                    peripheral_volume=TRUE_V2_L,
                    peripheral_volume_unit="L",
                    time=float(minute),
                    time_unit="min",
                    concentration_unit="mg/L",
                ),
                6,
            ),
        )
        for idx, (minute, mult) in enumerate(zip(quantified_minutes, _NOISY_MULTIPLIERS, strict=True), start=1)
    )
    return ConcentrationProfile(
        subject_id="SUBJ-NOISY-2C",
        analyte_id="PW-0002-PARENT",
        dose_event_id="DOSE-NOISY-2C",
        route=Route.IV_BOLUS,
        dose_mg=DOSE_MG,
        points=points,
        n_blq_omitted=1,
        n_above_uloq_omitted=0,
        n_predose_omitted=1,
    )


def main() -> None:
    """Write the synthetic two-compartment IV bolus study directory."""

    write_study_directory(build_two_compartment_study(), TWO_COMP_EXAMPLE_DIR)
    print(f"Wrote {TWO_COMP_EXAMPLE_DIR}")


if __name__ == "__main__":
    main()
