"""Concentration profiles for the locked single-dose IV bolus design.

NCA and the one-compartment fit both read this view. Building it does not
calculate parameters. Sparse sampling, multiple dosing, and non-bolus routes
are rejected here so they are not passed to a later calculator.
"""

from collections import defaultdict
from dataclasses import dataclass

from pharmacoweave.core.entities import (
    Assay,
    Biosample,
    DoseEvent,
    Observation,
    ObservationStatus,
    ObservationType,
    Route,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import convert
from pharmacoweave.errors import UnitError, UnsupportedDesign


@dataclass(frozen=True)
class TimedConcentration:
    """One sample on a concentration profile, already on the hour scale."""

    time_h: float
    status: ObservationStatus
    value: float | None
    value_unit: str | None
    observation_id: str
    is_predose: bool


@dataclass(frozen=True)
class ConcentrationProfile:
    """Serial concentrations for one subject, analyte, and dose event."""

    subject_id: str
    analyte_id: str
    dose_event_id: str
    route: Route
    dose_value: float
    dose_unit: str
    concentration_unit: str | None
    lloq: float | None
    points: tuple[TimedConcentration, ...]
    unit_error: str | None = None


def load_iv_bolus_profiles(study: PharmacologyStudy) -> tuple[ConcentrationProfile, ...]:
    """Return IV bolus profiles, or reject a design outside that lock.

    Rich serial sampling means at least one profile has two distinct post-dose
    times. A study in which every subject contributes a single post-dose time
    is sparse and is not analyzed by this calculator.
    """

    _require_single_dose_iv_bolus(study)
    profiles = _profiles(study)
    if not profiles:
        raise UnsupportedDesign("No concentration observations were found.")
    if not any(_postdose_time_count(profile) >= 2 for profile in profiles):
        raise UnsupportedDesign(
            "Sparse sampling is outside the locked rich_serial design. "
            "Every profile has fewer than two distinct post-dose times."
        )
    return tuple(sorted(profiles, key=lambda item: (item.subject_id, item.analyte_id, item.dose_event_id)))


def _require_single_dose_iv_bolus(study: PharmacologyStudy) -> None:
    if not study.dose_events:
        raise UnsupportedDesign("The locked IV bolus analysis needs at least one dose event.")
    by_subject: dict[str, list[DoseEvent]] = defaultdict(list)
    for dose in study.dose_events:
        by_subject[dose.subject_id].append(dose)
    for subject_id, doses in by_subject.items():
        if len(doses) != 1:
            raise UnsupportedDesign(
                f"{subject_id} has {len(doses)} dose events. "
                "Multiple dosing is outside the locked single_dose design."
            )
        route = doses[0].route
        if route is not Route.IV_BOLUS:
            raise UnsupportedDesign(
                f"{subject_id} was dosed by {route.value}. "
                "Only single-dose iv_bolus rich sampling is implemented."
            )


def _profiles(study: PharmacologyStudy) -> list[ConcentrationProfile]:
    samples = {sample.sample_id: sample for sample in study.biosamples}
    doses = {dose.dose_event_id: dose for dose in study.dose_events}
    assays = {assay.assay_id: assay for assay in study.assays}
    grouped: dict[tuple[str, str, str], list[tuple[Biosample, Observation]]] = defaultdict(list)
    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None or observation.analyte_id is None:
            raise UnsupportedDesign(
                f"Observation {observation.observation_id} is not linked to a sample and analyte."
            )
        if sample.dose_event_id not in doses:
            raise UnsupportedDesign(f"Sample {sample.sample_id} points at a missing dose event.")
        grouped[(sample.subject_id, observation.analyte_id, sample.dose_event_id)].append((sample, observation))

    profiles: list[ConcentrationProfile] = []
    for (subject_id, analyte_id, dose_event_id), rows in grouped.items():
        dose = doses[dose_event_id]
        profiles.append(
            _one_profile(
                subject_id=subject_id,
                analyte_id=analyte_id,
                dose_event_id=dose_event_id,
                dose=dose,
                rows=rows,
                assays=assays,
            )
        )
    return profiles


def _one_profile(
    *,
    subject_id: str,
    analyte_id: str,
    dose_event_id: str,
    dose: DoseEvent,
    rows: list[tuple[Biosample, Observation]],
    assays: dict[str, Assay],
) -> ConcentrationProfile:
    ordered = sorted(rows, key=lambda pair: _time_h(pair[0]))
    target_unit, unit_error = _common_unit(ordered)
    lloq = _lloq(ordered, assays, target_unit)
    if target_unit is None:
        target_unit = _assay_limit_unit(ordered, assays)
    points: list[TimedConcentration] = []
    for sample, observation in ordered:
        value = None
        value_unit = None
        if (
            observation.status is ObservationStatus.QUANTIFIED
            and unit_error is None
            and observation.numeric_value is not None
            and observation.unit is not None
            and target_unit is not None
        ):
            value = convert(observation.numeric_value, observation.unit, target_unit)
            value_unit = target_unit
        points.append(
            TimedConcentration(
                time_h=_time_h(sample),
                status=observation.status,
                value=value,
                value_unit=value_unit,
                observation_id=observation.observation_id,
                is_predose=sample.is_predose,
            )
        )
    return ConcentrationProfile(
        subject_id=subject_id,
        analyte_id=analyte_id,
        dose_event_id=dose_event_id,
        route=dose.route,
        dose_value=dose.dose_value,
        dose_unit=dose.dose_unit,
        concentration_unit=target_unit,
        lloq=lloq,
        points=tuple(points),
        unit_error=unit_error,
    )


def _time_h(sample: Biosample) -> float:
    return convert(sample.time_since_dose, sample.time_since_dose_unit, "h")


def _common_unit(rows: list[tuple[Biosample, Observation]]) -> tuple[str | None, str | None]:
    quantified = [observation for _sample, observation in rows if observation.status is ObservationStatus.QUANTIFIED]
    if not quantified:
        return None, None
    target = quantified[0].unit
    if target is None:
        return None, "A quantified concentration is missing its unit."
    for observation in quantified:
        if observation.unit is None or observation.numeric_value is None:
            return None, "A quantified concentration is missing its value or unit."
        try:
            convert(observation.numeric_value, observation.unit, target)
        except UnitError as exc:
            return None, f"Quantified concentrations do not share a convertible unit ({exc})."
    return target, None


def _lloq(
    rows: list[tuple[Biosample, Observation]],
    assays: dict[str, Assay],
    target_unit: str | None,
) -> float | None:
    if target_unit is None:
        return None
    for _sample, observation in rows:
        if observation.assay_id is None:
            continue
        assay = assays.get(observation.assay_id)
        if assay is None or assay.lloq is None or assay.limit_unit is None:
            continue
        try:
            return convert(assay.lloq, assay.limit_unit, target_unit)
        except UnitError:
            return None
    return None


def _assay_limit_unit(rows: list[tuple[Biosample, Observation]], assays: dict[str, Assay]) -> str | None:
    for _sample, observation in rows:
        if observation.assay_id is None:
            continue
        assay = assays.get(observation.assay_id)
        if assay is not None and assay.limit_unit:
            return assay.limit_unit
    return None


def _postdose_time_count(profile: ConcentrationProfile) -> int:
    times = {round(point.time_h, 9) for point in profile.points if point.time_h >= 0}
    return len(times)
