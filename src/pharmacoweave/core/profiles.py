"""Post-dose concentration series for one subject and one analyte.

BLQ and above-ULOQ rows are counted and left out of the numeric series.
They are not replaced with zero. Predose rows are not part of the series.
"""

from dataclasses import dataclass

from pharmacoweave.core.entities import Biosample, Observation, ObservationStatus, ObservationType, Route
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import CONC_MASS, DRUG_MASS, convert, parse_unit
from pharmacoweave.errors import AnalysisError


@dataclass(frozen=True)
class QuantifiedPoint:
    """One quantified concentration at a time since dose."""

    observation_id: str
    time_h: float
    concentration_mg_l: float


@dataclass(frozen=True)
class ConcentrationProfile:
    """The quantified curve used by single-dose IV bolus NCA and fitting."""

    subject_id: str
    analyte_id: str
    dose_event_id: str
    route: Route
    dose_mg: float
    points: tuple[QuantifiedPoint, ...]
    n_blq_omitted: int
    n_above_uloq_omitted: int
    n_predose_omitted: int


def iv_bolus_profiles(study: PharmacologyStudy) -> list[ConcentrationProfile]:
    """Build one single-dose IV bolus profile per subject and analyte.

    A subject with more than one dose event is rejected. This calculator does
    not superpose multiple doses.
    """

    doses_by_subject: dict[str, list[str]] = {}
    for dose in study.dose_events:
        doses_by_subject.setdefault(dose.subject_id, []).append(dose.dose_event_id)
    for subject_id, dose_ids in doses_by_subject.items():
        if len(dose_ids) != 1:
            raise AnalysisError(
                f"Subject {subject_id} has {len(dose_ids)} dose events. "
                "The implemented calculator accepts one single dose per subject."
            )
    doses = {dose.dose_event_id: dose for dose in study.dose_events}
    samples = {sample.sample_id: sample for sample in study.biosamples}
    analytes = {analyte.analyte_id: analyte for analyte in study.analytes}
    grouped: dict[tuple[str, str, str], list[tuple[Biosample, Observation]]] = {}
    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION:
            continue
        if observation.analyte_id is None:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None:
            raise AnalysisError(f"Observation {observation.observation_id} has no biosample.")
        dose_event = doses.get(sample.dose_event_id)
        if dose_event is None:
            raise AnalysisError(f"Sample {sample.sample_id} has no dose event.")
        grouped.setdefault((sample.subject_id, observation.analyte_id, dose_event.dose_event_id), []).append(
            (sample, observation)
        )
    profiles: list[ConcentrationProfile] = []
    for (subject_id, analyte_id, dose_event_id), rows in sorted(grouped.items()):
        dose = doses[dose_event_id]
        if dose.route is not Route.IV_BOLUS:
            raise AnalysisError(
                f"Dose {dose.dose_event_id} uses route {dose.route.value}. This calculator accepts iv_bolus only."
            )
        if parse_unit(dose.dose_unit).dimensions != DRUG_MASS:
            raise AnalysisError(f"Dose {dose.dose_event_id} must be a drug mass, such as mg.")
        analyte = analytes.get(analyte_id)
        if analyte is not None and analyte.compound_id is not None and analyte.compound_id != dose.compound_id:
            raise AnalysisError(
                f"Analyte {analyte_id} is linked to {analyte.compound_id}, not dose compound {dose.compound_id}."
            )
        points: list[QuantifiedPoint] = []
        n_blq = 0
        n_above = 0
        n_predose = 0
        for sample, observation in rows:
            time_h = convert(sample.time_since_dose, sample.time_since_dose_unit, "h")
            if sample.is_predose or time_h < 0:
                n_predose += 1
                continue
            if observation.status is ObservationStatus.BLQ:
                n_blq += 1
                continue
            if observation.status is ObservationStatus.ABOVE_ULOQ:
                n_above += 1
                continue
            if observation.status is not ObservationStatus.QUANTIFIED:
                continue
            if observation.unit is None or observation.numeric_value is None:
                raise AnalysisError(f"Quantified observation {observation.observation_id} has no value.")
            if parse_unit(observation.unit).dimensions != CONC_MASS:
                raise AnalysisError(
                    f"Observation {observation.observation_id} must be a mass concentration. "
                    "Molar units need an explicit molecular-weight conversion."
                )
            points.append(
                QuantifiedPoint(
                    observation_id=observation.observation_id,
                    time_h=time_h,
                    concentration_mg_l=convert(observation.numeric_value, observation.unit, "mg/L"),
                )
            )
        points.sort(key=lambda item: (item.time_h, item.observation_id))
        for earlier, later in zip(points, points[1:], strict=False):
            if abs(later.time_h - earlier.time_h) <= 1e-9:
                raise AnalysisError(
                    f"Subject {subject_id} has two quantified concentrations at the same time. "
                    "Collapse replicates before NCA."
                )
        profiles.append(
            ConcentrationProfile(
                subject_id=subject_id,
                analyte_id=analyte_id,
                dose_event_id=dose_event_id,
                route=dose.route,
                dose_mg=convert(dose.dose_value, dose.dose_unit, "mg"),
                points=tuple(points),
                n_blq_omitted=n_blq,
                n_above_uloq_omitted=n_above,
                n_predose_omitted=n_predose,
            )
        )
    return profiles
