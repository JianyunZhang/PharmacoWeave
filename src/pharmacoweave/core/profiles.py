"""Post-dose concentration series and auditable replicate collapsing.

BLQ and above-ULOQ rows are counted and left out of the numeric series.
They are not replaced with zero. Predose rows are not part of the series.
"""

import math
from collections.abc import Collection
from dataclasses import dataclass
from datetime import datetime
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import (
    Biosample,
    Compound,
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
from pharmacoweave.core.units import (
    AMOUNT,
    CONC_AMOUNT,
    CONC_MASS,
    DRUG_MASS,
    convert,
    convert_with_molecular_weight,
    parse_unit,
)
from pharmacoweave.errors import AnalysisError, UnitError

ReplicateMethod = Literal["arithmetic_mean", "geometric_mean"]
ReplicatePolicy = Literal["error", "arithmetic_mean", "geometric_mean"]

EXTRAVASCULAR_ROUTES: tuple[Route, ...] = (
    Route.ORAL,
    Route.SUBCUTANEOUS,
    Route.INTRAMUSCULAR,
)
SUPPORTED_SINGLE_DOSE_ROUTES: tuple[Route, ...] = (
    Route.IV_BOLUS,
    Route.IV_INFUSION,
    *EXTRAVASCULAR_ROUTES,
)


@dataclass(frozen=True)
class QuantifiedPoint:
    """One quantified concentration at a time since dose."""

    observation_id: str
    time_h: float
    concentration_mg_l: float


@dataclass(frozen=True)
class ConcentrationProfile:
    """The quantified curve used by single-dose NCA and compartmental fitting."""

    subject_id: str
    analyte_id: str
    dose_event_id: str
    route: Route
    dose_mg: float
    points: tuple[QuantifiedPoint, ...]
    n_blq_omitted: int
    n_above_uloq_omitted: int
    n_predose_omitted: int
    infusion_duration_h: float | None = None
    tau_h: float | None = None
    n_doses: int = 1
    interval_regimen: Literal["single_dose", "repeated_equal", "steady_state"] = "single_dose"
    reference_cmax_mg_l: float | None = None


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
    analytes = {analyte.analyte_id: analyte for analyte in study.analytes}
    compounds = {compound.compound_id: compound for compound in study.compounds}
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
        analyte = analytes.get(analyte_id)
        compound = compounds.get(analyte.compound_id) if analyte and analyte.compound_id else None
        first = group[0]
        if first.unit is None:
            raise AnalysisError(f"Quantified observation {first.observation_id} has no unit.")
        target_unit = first.unit
        values: list[float] = []
        for obs in group:
            if obs.numeric_value is None or obs.unit is None:
                raise AnalysisError(f"Quantified observation {obs.observation_id} has no value.")
            try:
                val = _convert_with_optional_compound(obs.numeric_value, obs.unit, target_unit, compound)
            except UnitError as exc:
                raise AnalysisError(str(exc)) from exc
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


def single_dose_profiles(
    study: PharmacologyStudy,
    *,
    allowed_routes: Collection[Route] | None = None,
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build one single-dose concentration profile per subject and analyte.

    A subject with more than one dose event is rejected. Duplicate same-time
    quantified concentrations raise an error unless ``replicate_policy`` is
    explicitly set to ``'arithmetic_mean'`` or ``'geometric_mean'``. Molar
    doses and concentrations are converted to ``mg`` and ``mg/L`` only when
    the linked ``Compound`` declares an explicit ``molecular_weight``.
    """

    if replicate_policy not in {"error", "arithmetic_mean", "geometric_mean"}:
        raise AnalysisError("replicate_policy must be 'error', 'arithmetic_mean', or 'geometric_mean'.")
    if replicate_policy != "error":
        study, _ = collapse_replicates(study, method=replicate_policy)

    permitted = frozenset(SUPPORTED_SINGLE_DOSE_ROUTES if allowed_routes is None else allowed_routes)
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
    compounds = {compound.compound_id: compound for compound in study.compounds}

    grouped: dict[tuple[str, str, str], list[tuple[Biosample, Observation]]] = {}
    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION:
            continue
        if observation.analyte_id is None:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None:
            raise AnalysisError(f"Observation {observation.observation_id} has no biosample.")
        if sample.matrix == "urine":
            continue
        dose_event = doses.get(sample.dose_event_id)
        if dose_event is None:
            raise AnalysisError(f"Sample {sample.sample_id} has no dose event.")
        grouped.setdefault((sample.subject_id, observation.analyte_id, dose_event.dose_event_id), []).append(
            (sample, observation)
        )

    profiles: list[ConcentrationProfile] = []
    for (subject_id, analyte_id, dose_event_id), rows in sorted(grouped.items()):
        dose = doses[dose_event_id]
        if dose.route not in permitted:
            allowed_label = ", ".join(sorted(route.value for route in permitted))
            raise AnalysisError(
                f"Dose {dose.dose_event_id} uses route {dose.route.value}. "
                f"This calculator accepts {allowed_label} only."
            )
        analyte = analytes.get(analyte_id)
        if analyte is not None and analyte.compound_id is not None and analyte.compound_id != dose.compound_id:
            raise AnalysisError(
                f"Analyte {analyte_id} is linked to {analyte.compound_id}, not dose compound {dose.compound_id}."
            )
        dose_compound = compounds.get(dose.compound_id)
        analyte_compound = (
            compounds.get(analyte.compound_id) if analyte is not None and analyte.compound_id is not None else None
        )
        dose_mg = _dose_to_mg(dose.dose_event_id, dose.dose_value, dose.dose_unit, dose_compound)
        infusion_duration_h: float | None = None
        if dose.route is Route.IV_INFUSION:
            if dose.infusion_duration is None or dose.infusion_duration_unit is None:
                raise AnalysisError(f"IV infusion dose {dose.dose_event_id} requires infusion_duration.")
            infusion_duration_h = convert(dose.infusion_duration, dose.infusion_duration_unit, "h")
            if infusion_duration_h <= 0:
                raise AnalysisError(f"IV infusion dose {dose.dose_event_id} requires a positive infusion_duration.")

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
            conc_mg_l = _concentration_to_mg_l(
                observation.observation_id,
                observation.numeric_value,
                observation.unit,
                analyte_compound or dose_compound,
            )
            points.append(
                QuantifiedPoint(
                    observation_id=observation.observation_id,
                    time_h=time_h,
                    concentration_mg_l=conc_mg_l,
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
                dose_mg=dose_mg,
                points=tuple(points),
                n_blq_omitted=n_blq,
                n_above_uloq_omitted=n_above,
                n_predose_omitted=n_predose,
                infusion_duration_h=infusion_duration_h,
            )
        )
    return profiles


def iv_bolus_profiles(
    study: PharmacologyStudy,
    *,
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build one single-dose IV bolus profile per subject and analyte."""

    return single_dose_profiles(
        study,
        allowed_routes=(Route.IV_BOLUS,),
        replicate_policy=replicate_policy,
    )


def iv_infusion_profiles(
    study: PharmacologyStudy,
    *,
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build one single-dose constant-rate IV infusion profile per subject and analyte."""

    return single_dose_profiles(
        study,
        allowed_routes=(Route.IV_INFUSION,),
        replicate_policy=replicate_policy,
    )


def extravascular_profiles(
    study: PharmacologyStudy,
    *,
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build one single-dose extravascular (oral, subcutaneous, intramuscular) profile per subject and analyte."""

    return single_dose_profiles(
        study,
        allowed_routes=EXTRAVASCULAR_ROUTES,
        replicate_policy=replicate_policy,
    )


def equal_interval_profiles(
    study: PharmacologyStudy,
    *,
    allowed_routes: Collection[Route] | None = None,
    regimen: Literal["repeated_equal", "steady_state"],
    replicate_policy: ReplicatePolicy = "error",
) -> list[ConcentrationProfile]:
    """Build the analyzed interdose interval for equal doses and equal spacing.

    ``repeated_equal`` uses the last interval that ends at the next dose.
    ``steady_state`` uses the interval that starts at the last dose and lasts
    one spacing. Samples are placed by absolute clock time. A sample exactly
    at the end of the interval is the end-of-interval concentration, not the
    concentration after the next dose.
    """

    if regimen not in {"repeated_equal", "steady_state"}:
        raise AnalysisError("regimen must be repeated_equal or steady_state.")
    if replicate_policy not in {"error", "arithmetic_mean", "geometric_mean"}:
        raise AnalysisError("replicate_policy must be 'error', 'arithmetic_mean', or 'geometric_mean'.")
    if replicate_policy != "error":
        study, _ = collapse_replicates(study, method=replicate_policy)

    permitted = frozenset(SUPPORTED_SINGLE_DOSE_ROUTES if allowed_routes is None else allowed_routes)
    by_subject: dict[str, list] = {}
    for dose in study.dose_events:
        by_subject.setdefault(dose.subject_id, []).append(dose)
    samples = {sample.sample_id: sample for sample in study.biosamples}
    analytes = {analyte.analyte_id: analyte for analyte in study.analytes}
    compounds = {compound.compound_id: compound for compound in study.compounds}

    profiles: list[ConcentrationProfile] = []
    grouped: dict[tuple[str, str], list[tuple[Biosample, Observation]]] = {}
    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION or observation.analyte_id is None:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None:
            raise AnalysisError(f"Observation {observation.observation_id} has no biosample.")
        if sample.matrix == "urine":
            continue
        grouped.setdefault((sample.subject_id, observation.analyte_id), []).append((sample, observation))

    for subject_id, subject_doses in sorted(by_subject.items()):
        ordered = sorted(subject_doses, key=lambda item: item.absolute_time)
        if len(ordered) < 2:
            raise AnalysisError(
                f"Subject {subject_id} has {len(ordered)} dose events. "
                "Equal-interval NCA requires at least two doses."
            )
        if any(dose.route not in permitted for dose in ordered):
            raise AnalysisError(f"Subject {subject_id} has a route outside this NCA design.")
        if len({dose.route for dose in ordered}) != 1 or len({dose.compound_id for dose in ordered}) != 1:
            raise AnalysisError(f"Subject {subject_id} mixes route or compound across doses.")
        compound = compounds.get(ordered[0].compound_id)
        dose_mg = [
            _dose_to_mg(dose.dose_event_id, dose.dose_value, dose.dose_unit, compound) for dose in ordered
        ]
        if any(abs(amount - dose_mg[0]) > 1e-9 * max(1.0, dose_mg[0]) for amount in dose_mg):
            raise AnalysisError(f"Subject {subject_id} doses are not equal. Unequal regimens are not calculated.")
        gaps = [
            (later.absolute_time - earlier.absolute_time).total_seconds() / 3600.0
            for earlier, later in zip(ordered, ordered[1:], strict=False)
        ]
        if any(gap <= 0 for gap in gaps) or any(abs(gap - gaps[0]) > 1e-6 for gap in gaps):
            raise AnalysisError(f"Subject {subject_id} dosing intervals are not equal and positive.")
        tau_h = gaps[0]
        infusions = [_infusion_hours(dose) for dose in ordered]
        if ordered[0].route is Route.IV_INFUSION:
            if any(value is None or value <= 0 for value in infusions):
                raise AnalysisError(f"Subject {subject_id} is missing a positive infusion duration.")
            assert infusions[0] is not None
            if any(value is None or abs(value - infusions[0]) > 1e-6 for value in infusions):
                raise AnalysisError(f"Subject {subject_id} infusion durations are not equal.")
            if infusions[0] >= tau_h:
                raise AnalysisError(
                    f"Subject {subject_id} infusion duration is not shorter than the dosing interval."
                )
        anchor = ordered[-1] if regimen == "steady_state" else ordered[-2]
        analyte_ids = {analyte_id for (sid, analyte_id) in grouped if sid == subject_id}
        if not analyte_ids:
            raise AnalysisError(f"Subject {subject_id} has no concentration observations.")
        for analyte_id in sorted(analyte_ids):
            analyte = analytes.get(analyte_id)
            if analyte is not None and analyte.compound_id not in {None, ordered[0].compound_id}:
                raise AnalysisError(
                    f"Analyte {analyte_id} is linked to {analyte.compound_id}, "
                    f"not dose compound {ordered[0].compound_id}."
                )
            analyte_compound = (
                compounds.get(analyte.compound_id) if analyte is not None and analyte.compound_id is not None else None
            )
            points, n_blq, n_above, n_predose = _interval_points(
                grouped[(subject_id, analyte_id)],
                anchor_time=anchor.absolute_time,
                tau_h=tau_h,
                compound=analyte_compound or compound,
            )
            reference = None
            if regimen == "repeated_equal" and len(ordered) >= 3:
                first_points, _, _, _ = _interval_points(
                    grouped[(subject_id, analyte_id)],
                    anchor_time=ordered[0].absolute_time,
                    tau_h=tau_h,
                    compound=analyte_compound or compound,
                )
                if first_points:
                    reference = max(point.concentration_mg_l for point in first_points)
            infusion_h = infusions[0] if ordered[0].route is Route.IV_INFUSION else None
            profiles.append(
                ConcentrationProfile(
                    subject_id=subject_id,
                    analyte_id=analyte_id,
                    dose_event_id=anchor.dose_event_id,
                    route=ordered[0].route,
                    dose_mg=dose_mg[0],
                    points=tuple(points),
                    n_blq_omitted=n_blq,
                    n_above_uloq_omitted=n_above,
                    n_predose_omitted=n_predose,
                    infusion_duration_h=infusion_h,
                    tau_h=tau_h,
                    n_doses=len(ordered),
                    interval_regimen=regimen,
                    reference_cmax_mg_l=reference,
                )
            )
    return profiles


def _infusion_hours(dose) -> float | None:
    if dose.infusion_duration is None or dose.infusion_duration_unit is None:
        return None
    return convert(dose.infusion_duration, dose.infusion_duration_unit, "h")


def _interval_points(
    rows: list[tuple[Biosample, Observation]],
    *,
    anchor_time: datetime,
    tau_h: float,
    compound: Compound | None,
) -> tuple[list[QuantifiedPoint], int, int, int]:
    points: list[QuantifiedPoint] = []
    n_blq = 0
    n_above = 0
    n_predose = 0
    for sample, observation in rows:
        time_h = (sample.collection_time_absolute - anchor_time).total_seconds() / 3600.0
        if time_h < -1e-6:
            n_predose += 1
            continue
        if time_h > tau_h + 1e-6:
            continue
        if abs(time_h) <= 1e-6:
            time_h = 0.0
        elif abs(time_h - tau_h) <= 1e-6:
            time_h = tau_h
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
        points.append(
            QuantifiedPoint(
                observation_id=observation.observation_id,
                time_h=time_h,
                concentration_mg_l=_concentration_to_mg_l(
                    observation.observation_id,
                    observation.numeric_value,
                    observation.unit,
                    compound,
                ),
            )
        )
    points.sort(key=lambda item: (item.time_h, item.observation_id))
    for earlier, later in zip(points, points[1:], strict=False):
        if abs(later.time_h - earlier.time_h) <= 1e-9:
            raise AnalysisError(
                f"Two quantified concentrations fall at {earlier.time_h:g} h in one dosing interval. "
                "Collapse replicates before NCA."
            )
    return points, n_blq, n_above, n_predose


def _convert_with_optional_compound(
    value: float,
    source: str,
    target: str,
    compound: Compound | None,
) -> float:
    return convert_with_molecular_weight(
        value,
        source,
        target,
        molecular_weight=None if compound is None else compound.molecular_weight,
        molecular_weight_unit=None if compound is None else compound.molecular_weight_unit,
    )


def _dose_to_mg(dose_event_id: str, value: float, unit: str, compound: Compound | None) -> float:
    dimensions = parse_unit(unit).dimensions
    if dimensions == DRUG_MASS:
        return convert(value, unit, "mg")
    if dimensions == AMOUNT:
        if compound is None or compound.molecular_weight is None or compound.molecular_weight_unit is None:
            raise AnalysisError(
                f"Dose {dose_event_id} uses molar unit {unit!r}. "
                "Converting molar doses requires Compound.molecular_weight and molecular_weight_unit."
            )
        try:
            return _convert_with_optional_compound(value, unit, "mg", compound)
        except UnitError as exc:
            raise AnalysisError(str(exc)) from exc
    raise AnalysisError(f"Dose {dose_event_id} must be a drug mass or molar amount, such as mg or umol.")


def _concentration_to_mg_l(
    observation_id: str,
    value: float,
    unit: str,
    compound: Compound | None,
) -> float:
    dimensions = parse_unit(unit).dimensions
    if dimensions == CONC_MASS:
        return convert(value, unit, "mg/L")
    if dimensions == CONC_AMOUNT:
        if compound is None or compound.molecular_weight is None or compound.molecular_weight_unit is None:
            raise AnalysisError(
                f"Observation {observation_id} uses molar concentration {unit!r}. "
                "Molar units need an explicit molecular-weight conversion on Compound."
            )
        try:
            return _convert_with_optional_compound(value, unit, "mg/L", compound)
        except UnitError as exc:
            raise AnalysisError(str(exc)) from exc
    raise AnalysisError(f"Observation {observation_id} must be a mass or molar concentration.")
