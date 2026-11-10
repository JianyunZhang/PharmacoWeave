"""Urine excretion amount and renal clearance.

Interval amount is urine concentration times recorded void volume. The plasma
trapezoidal rule is not applied to urine concentrations. A missing volume or a
BLQ urine concentration leaves the amount empty. Renal clearance is the summed
amount divided by a caller-supplied plasma AUC, and only when every interval
in that collection has an amount.
"""

import math
from pathlib import Path
from typing import Literal

import yaml
from pydantic import Field, ValidationError, model_validator

from pharmacoweave.core.entities import (
    Biosample,
    Compound,
    Observation,
    ObservationStatus,
    ObservationType,
    Record,
    UrineCollection,
)
from pharmacoweave.core.profiles import _concentration_to_mg_l
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import AUC_MASS, VOLUME, convert, parse_unit
from pharmacoweave.core.validation import validate_study
from pharmacoweave.errors import AnalysisError, IngestError, ValidationFailed
from pharmacoweave.version import SCHEMA_VERSION

_LIMITATIONS = (
    "Urine amount is concentration times recorded void volume.",
    "BLQ urine concentrations and missing void volumes stay empty and are not replaced with zero.",
    "Urine concentration is not a plasma concentration and is not integrated by the plasma trapezoidal method.",
    "Renal clearance is the summed urine amount divided by a caller-supplied plasma AUC for the same window.",
    "Research use only. The estimates are not a dosing recommendation.",
)


class UrineNcaConfig(Record):
    """Rules for urine excretion. This config has no plasma AUC method."""

    analysis: Literal["urine_excretion"]
    schema_version: Literal["0.1.0"] = "0.1.0"
    note: str = ""

    @model_validator(mode="after")
    def version_matches_package(self) -> "UrineNcaConfig":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        return self


class PlasmaAucInput(Record):
    """A plasma AUC supplied by the caller. It is not calculated from urine."""

    source_id: str = Field(min_length=1)
    subject_id: str = Field(min_length=1)
    analyte_id: str = Field(min_length=1)
    auc: float = Field(gt=0)
    auc_unit: str = Field(min_length=1)
    matrix: Literal["plasma", "serum", "whole_blood", "blood"]

    @model_validator(mode="after")
    def auc_is_plasma_exposure(self) -> "PlasmaAucInput":
        if parse_unit(self.auc_unit).dimensions != AUC_MASS:
            raise ValueError("auc_unit must be a mass-time/volume AUC unit such as mg*h/L.")
        return self


class UrineIntervalAmount(Record):
    """One collection interval. Empty amounts stay null."""

    sample_id: str
    observation_id: str
    collection_id: str | None = None
    interval_start_h: float | None = None
    interval_end_h: float | None = None
    urine_concentration_mg_l: float | None = None
    urine_volume_l: float | None = None
    amount_excreted_mg: float | None = None
    omitted_reason: str | None = None


class UrineExcretionResult(Record):
    """Subject-level urine excretion. Plasma AUC columns are not stored here."""

    subject_id: str
    analyte_id: str
    dose_event_id: str
    matrix: Literal["urine"] = "urine"
    intervals: list[UrineIntervalAmount] = Field(default_factory=list)
    amount_excreted_mg: float | None = None
    renal_clearance_l_h: float | None = None
    plasma_auc_source_id: str | None = None
    limitations: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def clearance_needs_a_complete_amount(self) -> "UrineExcretionResult":
        if self.renal_clearance_l_h is not None and self.amount_excreted_mg is None:
            raise ValueError("Renal clearance requires a complete excreted amount.")
        if self.matrix != "urine":
            raise ValueError("Urine excretion results stay on matrix urine.")
        return self


def load_urine_nca_config(path: str | Path) -> UrineNcaConfig:
    """Load a urine excretion config. A plasma NCA config is not accepted."""

    selected = Path(path)
    if not selected.is_file():
        raise IngestError(f"Urine NCA config {selected} does not exist.")
    try:
        payload = yaml.safe_load(selected.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise IngestError(f"{selected} is not valid YAML: {exc}") from exc
    if not isinstance(payload, dict):
        raise IngestError(f"{selected} must contain a mapping.")
    try:
        return UrineNcaConfig.model_validate(payload)
    except ValidationError as exc:
        raise IngestError(f"{selected} is not a urine excretion config: {exc}") from exc


def load_plasma_auc_inputs(path: str | Path) -> list[PlasmaAucInput]:
    """Load caller-supplied plasma AUC values. Urine matrix is rejected."""

    selected = Path(path)
    if not selected.is_file():
        raise IngestError(f"Plasma AUC file {selected} does not exist.")
    try:
        payload = yaml.safe_load(selected.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise IngestError(f"{selected} is not valid YAML: {exc}") from exc
    rows = payload if isinstance(payload, list) else [payload]
    if not isinstance(payload, list | dict):
        raise IngestError(f"{selected} must contain a mapping or a list of mappings.")
    loaded: list[PlasmaAucInput] = []
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            raise IngestError(f"{selected} item {index} must be a mapping.")
        try:
            loaded.append(PlasmaAucInput.model_validate(row))
        except ValidationError as exc:
            raise IngestError(f"{selected} item {index} is not a plasma AUC: {exc}") from exc
    return loaded


def run_urine_excretion(
    study: PharmacologyStudy,
    config: str | Path | UrineNcaConfig,
    *,
    plasma_aucs: list[PlasmaAucInput] | None = None,
) -> list[UrineExcretionResult]:
    """Calculate urine interval amounts and, when complete, renal clearance."""

    if not isinstance(config, UrineNcaConfig):
        load_urine_nca_config(config)
    report = validate_study(study)
    if not report.ok:
        raise ValidationFailed(report)
    auc_by_key = {(item.subject_id, item.analyte_id): item for item in plasma_aucs or []}
    grouped: dict[tuple[str, str, str], list[tuple[Biosample, Observation]]] = {}
    samples = {sample.sample_id: sample for sample in study.biosamples}
    for observation in study.observations:
        if observation.observation_type is not ObservationType.CONCENTRATION or observation.analyte_id is None:
            continue
        sample = samples.get(observation.sample_id)
        if sample is None or sample.matrix != "urine":
            continue
        grouped.setdefault((sample.subject_id, observation.analyte_id, sample.dose_event_id), []).append(
            (sample, observation)
        )
    compounds = {item.compound_id: item for item in study.compounds}
    analytes = {item.analyte_id: item for item in study.analytes}
    collections = {item.sample_id: item for item in study.urine_collections}
    results: list[UrineExcretionResult] = []
    for key in sorted(grouped):
        subject_id, analyte_id, dose_event_id = key
        analyte = analytes.get(analyte_id)
        compound = None
        if analyte is not None and analyte.compound_id is not None:
            compound = compounds.get(analyte.compound_id)
        intervals = _intervals(grouped[key], collections, compound)
        _reject_overlap(subject_id, intervals)
        amount, warnings = _complete_amount(intervals)
        auc = auc_by_key.get((subject_id, analyte_id))
        clearance = None
        source_id = None
        if amount is not None and auc is not None:
            auc_mg_h_l = convert(auc.auc, auc.auc_unit, "mg*h/L")
            clearance = amount / auc_mg_h_l
            source_id = auc.source_id
        elif amount is not None and auc is None:
            warnings.append("Renal clearance was withheld because no plasma AUC was supplied for this subject.")
        results.append(
            UrineExcretionResult(
                subject_id=subject_id,
                analyte_id=analyte_id,
                dose_event_id=dose_event_id,
                intervals=intervals,
                amount_excreted_mg=amount,
                renal_clearance_l_h=clearance,
                plasma_auc_source_id=source_id,
                limitations=list(_LIMITATIONS),
                warnings=warnings,
            )
        )
    return results


def _intervals(
    rows: list[tuple[Biosample, Observation]],
    collections: dict[str, UrineCollection],
    compound: Compound | None,
) -> list[UrineIntervalAmount]:
    quantified: dict[str, Observation] = {}
    chosen: dict[str, Observation] = {}
    for sample, observation in rows:
        previous = chosen.get(sample.sample_id)
        if previous is not None and previous.status is ObservationStatus.QUANTIFIED:
            if observation.status is ObservationStatus.QUANTIFIED:
                raise AnalysisError(
                    f"Sample {sample.sample_id} has two quantified urine concentrations. Collapse replicates first."
                )
            continue
        if observation.status is ObservationStatus.QUANTIFIED or previous is None:
            chosen[sample.sample_id] = observation
        if observation.status is ObservationStatus.QUANTIFIED:
            quantified[sample.sample_id] = observation
    built: list[UrineIntervalAmount] = []
    for sample, _observation in rows:
        if sample.sample_id in {item.sample_id for item in built}:
            continue
        observation = quantified.get(sample.sample_id, chosen[sample.sample_id])
        collection = collections.get(sample.sample_id)
        built.append(_one_interval(sample, observation, collection, compound))
    built.sort(key=lambda item: (item.interval_end_h is None, item.interval_end_h or 0.0, item.sample_id))
    return built


def _one_interval(
    sample: Biosample,
    observation: Observation,
    collection: UrineCollection | None,
    compound: Compound | None,
) -> UrineIntervalAmount:
    end_h = convert(sample.time_since_dose, sample.time_since_dose_unit, "h")
    start_h = None
    volume_l = None
    if collection is not None:
        start_h = convert(collection.interval_start_since_dose, collection.interval_start_unit, "h")
        if collection.volume is not None and collection.volume_unit is not None:
            if parse_unit(collection.volume_unit).dimensions != VOLUME:
                raise AnalysisError(f"Urine collection {collection.collection_id} volume unit is not a volume.")
            volume_l = convert(collection.volume, collection.volume_unit, "L")
    if observation.status is ObservationStatus.BLQ:
        return UrineIntervalAmount(
            sample_id=sample.sample_id,
            observation_id=observation.observation_id,
            collection_id=None if collection is None else collection.collection_id,
            interval_start_h=start_h,
            interval_end_h=end_h,
            omitted_reason="BLQ urine concentration was not replaced with zero.",
        )
    quantified = observation.status is ObservationStatus.QUANTIFIED
    if not quantified or observation.numeric_value is None or observation.unit is None:
        return UrineIntervalAmount(
            sample_id=sample.sample_id,
            observation_id=observation.observation_id,
            collection_id=None if collection is None else collection.collection_id,
            interval_start_h=start_h,
            interval_end_h=end_h,
            omitted_reason="The urine observation has no quantified concentration.",
        )
    concentration = _concentration_to_mg_l(
        observation.observation_id,
        observation.numeric_value,
        observation.unit,
        compound,
    )
    if collection is None or start_h is None:
        return UrineIntervalAmount(
            sample_id=sample.sample_id,
            observation_id=observation.observation_id,
            interval_end_h=end_h,
            urine_concentration_mg_l=concentration,
            omitted_reason="The urine collection interval was not recorded.",
        )
    if volume_l is None:
        return UrineIntervalAmount(
            sample_id=sample.sample_id,
            observation_id=observation.observation_id,
            collection_id=collection.collection_id,
            interval_start_h=start_h,
            interval_end_h=end_h,
            urine_concentration_mg_l=concentration,
            omitted_reason="Urine volume is missing and was not replaced with zero.",
        )
    return UrineIntervalAmount(
        sample_id=sample.sample_id,
        observation_id=observation.observation_id,
        collection_id=collection.collection_id,
        interval_start_h=start_h,
        interval_end_h=end_h,
        urine_concentration_mg_l=concentration,
        urine_volume_l=volume_l,
        amount_excreted_mg=concentration * volume_l,
    )


def _reject_overlap(subject_id: str, intervals: list[UrineIntervalAmount]) -> None:
    timed = [item for item in intervals if item.interval_start_h is not None and item.interval_end_h is not None]
    ordered = sorted(timed, key=lambda item: (item.interval_start_h or 0.0, item.interval_end_h or 0.0))
    for earlier, later in zip(ordered, ordered[1:], strict=False):
        if (later.interval_start_h or 0.0) < (earlier.interval_end_h or 0.0) - 1e-9:
            raise AnalysisError(f"Subject {subject_id} has overlapping urine collection intervals.")


def _complete_amount(intervals: list[UrineIntervalAmount]) -> tuple[float | None, list[str]]:
    warnings: list[str] = []
    if not intervals or any(item.amount_excreted_mg is None for item in intervals):
        warnings.append(
            "Total excretion and renal clearance were withheld because at least one interval has no amount."
        )
        return None, warnings
    ordered = sorted(intervals, key=lambda item: item.interval_start_h or 0.0)
    for earlier, later in zip(ordered, ordered[1:], strict=False):
        gap = (later.interval_start_h or 0.0) - (earlier.interval_end_h or 0.0)
        if gap > 1e-9:
            warnings.append(
                "Collection intervals are not contiguous, so total excretion and renal clearance were withheld."
            )
            return None, warnings
    total = math.fsum(item.amount_excreted_mg or 0.0 for item in ordered)
    return total, warnings
