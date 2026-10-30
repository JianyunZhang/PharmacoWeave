"""Post-dose concentration series and auditable replicate collapsing.

BLQ and above-ULOQ rows are counted and left out of the numeric series.
They are not replaced with zero. Predose rows are not part of the series.
"""

import math
from dataclasses import dataclass
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import (
    Biosample,
    FiniteFloat,
    Identifier,
    Observation,
    ObservationStatus,
    ObservationType,
    PositiveFloat,
    Record,
    Route,
)
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import CONC_MASS, DRUG_MASS, convert, parse_unit
from pharmacoweave.errors import AnalysisError

ReplicateMethod = Literal["arithmetic_mean", "geometric_mean"]
ReplicatePolicy = Literal["error", "arithmetic_mean", "geometric_mean"]


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


class ReplicateSummary(Record):
    """Audit record for collapsing quantified replicates at one post-dose time."""

    collapsed_observation_id: Identifier
    subject_id: Identifier
    analyte_id: Identifier
    dose_event_id: Identifier
    sample_id: Identifier
    time_h: FiniteFloat
    method: ReplicateMethod
    n_replicates: int = Field(ge=2)
    source_observation_ids: list[Identifier] = Field(min_length=2)
    collapsed_value: PositiveFloat
    unit: str
    cv_percent: float = Field(ge=0)
    n_censored_unchanged: int = Field(default=0, ge=0)


def collapse_replicates(
    study: PharmacologyStudy,
    *,
    method: ReplicateMethod = "arithmetic_mean",
) -> tuple[PharmacologyStudy, list[ReplicateSummary]]:
    """Collapse same-time post-dose quantified concentration replicates.

    BLQ, above-ULOQ, missing, invalid, and predose observations are never
    averaged or replaced with zero; they remain unchanged in the returned study.
    """

    if method not in {"arithmetic_mean", "geometric_mean"}:
        raise AnalysisError("method must be 'arithmetic_mean' or 'geometric_mean'.")

    samples = {sample.sample_id: sample for sample in study.biosamples}
    doses = {dose.dose_event_id: dose for dose in study.dose_events}
    existing_ids = {obs.observation_id for obs in study.observations}

    quantified_groups: dict[tuple[str, str, str, int], list[Observation]] = {}
    censored_counts: dict[tuple[str, str, str, int], int] = {}

    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION:
            continue
        if observation.analyte_id is None:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None:
            raise AnalysisError(f"Observation {observation.observation_id} has no biosample.")
        if sample.dose_event_id not in doses:
            raise AnalysisError(f"Sample {sample.sample_id} has no dose event.")
        time_h = convert(sample.time_since_dose, sample.time_since_dose_unit, "h")
        if sample.is_predose or time_h < 0:
            continue
        key = (sample.subject_id, observation.analyte_id, sample.dose_event_id, round(time_h * 1_000_000_000))
        if observation.status is ObservationStatus.QUANTIFIED:
            quantified_groups.setdefault(key, []).append(observation)
        elif observation.status in {ObservationStatus.BLQ, ObservationStatus.ABOVE_ULOQ}:
            censored_counts[key] = censored_counts.get(key, 0) + 1

    replacements: dict[str, Observation] = {}
    dropped_ids: set[str] = set()
    summaries: list[ReplicateSummary] = []

    for key in sorted(quantified_groups):
        group = sorted(quantified_groups[key], key=lambda item: item.observation_id)
        if len(group) < 2:
            continue
        subject_id, analyte_id, dose_event_id, _ = key
        forms = {obs.analyte_form for obs in group}
        if len(forms) > 1:
            raise AnalysisError(
                f"Subject {subject_id} has same-time observations with mixed analyte_form values; "
                "refusing to collapse distinct analyte forms."
            )
        first = group[0]
        if first.unit is None:
            raise AnalysisError(f"Quantified observation {first.observation_id} has no unit.")
        target_unit = first.unit
        values: list[float] = []
        for obs in group:
            if obs.numeric_value is None or obs.unit is None:
                raise AnalysisError(f"Quantified observation {obs.observation_id} has no value.")
            val = convert(obs.numeric_value, obs.unit, target_unit)
            if val <= 0:
                raise AnalysisError(f"Quantified concentration {obs.observation_id} must be positive.")
            values.append(val)

        k = len(values)
        arith_mean = sum(values) / k
        variance = sum((val - arith_mean) ** 2 for val in values) / (k - 1)
        sd = math.sqrt(max(0.0, variance))
        cv_percent = (100.0 * sd / arith_mean) if arith_mean > 0 else 0.0
        collapsed_val = (
            arith_mean
            if method == "arithmetic_mean"
            else math.exp(sum(math.log(val) for val in values) / k)
        )

        candidate_id = f"{first.observation_id}.collapsed"
        suffix = 2
        while candidate_id in existing_ids or len(candidate_id) > 128:
            base = first.observation_id[:96]
            candidate_id = f"{base}.rep{suffix}"
            suffix += 1
        existing_ids.add(candidate_id)

        source_ids = [obs.observation_id for obs in group]
        flags = [flag for flag in first.qc_flags if not flag.startswith("replicates_collapsed:")]
        flags.append(f"replicates_collapsed:{method}:n={k}")
        collapsed_obs = first.model_copy(
            update={
                "observation_id": candidate_id,
                "numeric_value": collapsed_val,
                "unit": target_unit,
                "raw_value": None,
                "reported_value": f"{method}({','.join(source_ids)})",
                "replicate_id": "collapsed",
                "qc_flags": flags,
            }
        )
        replacements[first.observation_id] = collapsed_obs
        for extra in group[1:]:
            dropped_ids.add(extra.observation_id)

        sample = samples[first.sample_id]
        time_h = convert(sample.time_since_dose, sample.time_since_dose_unit, "h")
        summaries.append(
            ReplicateSummary(
                collapsed_observation_id=candidate_id,
                subject_id=subject_id,
                analyte_id=analyte_id,
                dose_event_id=dose_event_id,
                sample_id=first.sample_id,
                time_h=time_h,
                method=method,
                n_replicates=k,
                source_observation_ids=source_ids,
                collapsed_value=collapsed_val,
                unit=target_unit,
                cv_percent=cv_percent,
                n_censored_unchanged=censored_counts.get(key, 0),
            )
        )

    if not summaries:
        return study, []

    updated_observations: list[Observation] = []
    for obs in study.observations:
        if obs.observation_id in dropped_ids:
            continue
        if obs.observation_id in replacements:
            updated_observations.append(replacements[obs.observation_id])
        else:
            updated_observations.append(obs)

    collapsed_study = study.model_copy(update={"observations": updated_observations})
    collapsed_study.attach_source_files(list(study.source_files))
    return collapsed_study, summaries


def iv_bolus_profiles(
    study: PharmacologyStudy,
    *,
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build one single-dose IV bolus profile per subject and analyte.

    A subject with more than one dose event is rejected. This calculator does
    not superpose multiple doses. Duplicate same-time quantified concentrations
    raise an error unless ``replicate_policy`` is explicitly set to
    ``'arithmetic_mean'`` or ``'geometric_mean'``.
    """

    if replicate_policy not in {"error", "arithmetic_mean", "geometric_mean"}:
        raise AnalysisError("replicate_policy must be 'error', 'arithmetic_mean', or 'geometric_mean'.")
    if replicate_policy != "error":
        study, _ = collapse_replicates(study, method=replicate_policy)

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
