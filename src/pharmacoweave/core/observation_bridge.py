"""Bridge study observation tables to typed EffectObservation and ToxEndpoint records.

Observations recorded in the study table with observation_type='effect' or
'toxicity' are converted into explicit downstream domain records with validated
dimensions, biosample/dose event relationships, and mandatory research-use
limitations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from pharmacoweave.core.study import PharmacologyStudy
    from pharmacoweave.pd_models.records import (
        EffectDimension,
        EffectObservation,
    )
    from pharmacoweave.toxicology.records import ToxEndpoint, ToxKind

from pharmacoweave.core.entities import (
    Assay,
    Biosample,
    DoseEvent,
    ObservationType,
)
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION, TOX_RESEARCH_LIMITATION
from pharmacoweave.core.units import (
    BODY_MASS,
    CONC_AMOUNT,
    CONC_MASS,
    DIMENSIONLESS,
    RATE,
    TIME,
    parse_unit,
)
from pharmacoweave.errors import IngestError, UnitError


def infer_effect_dimension(unit: str | None) -> EffectDimension:
    """Infer EffectDimension from a unit string."""
    from pharmacoweave.pd_models.records import EffectDimension

    if unit is None:
        return EffectDimension.FRACTION
    parsed = parse_unit(unit)
    if parsed.dimensions == DIMENSIONLESS:
        if parsed.factor == 0.01 or unit in {"%", "percent"}:
            return EffectDimension.PERCENT
        return EffectDimension.FRACTION
    if parsed.dimensions in {CONC_MASS, CONC_AMOUNT}:
        return EffectDimension.CONCENTRATION
    if parsed.dimensions == RATE:
        return EffectDimension.RATE
    if parsed.dimensions == TIME:
        return EffectDimension.TIME
    raise UnitError(f"Unit {unit!r} does not correspond to a supported effect dimension.")


def infer_tox_kind(matrix: str | None, unit: str | None) -> ToxKind:
    """Infer ToxKind from biological matrix and unit."""
    from pharmacoweave.toxicology.records import ToxKind

    if unit is not None:
        try:
            parsed = parse_unit(unit)
            if parsed.dimensions == BODY_MASS or (
                unit in {"percent", "fraction", "ratio", "%"}
                and matrix in {"whole_body", "body_weight"}
            ):
                return ToxKind.BODY_WEIGHT_CHANGE
            if parsed.dimensions in {CONC_MASS, CONC_AMOUNT, RATE}:
                return ToxKind.SERUM_CHEMISTRY
        except UnitError:
            pass
    if matrix is not None and matrix.lower() in {"whole_body", "body_weight"}:
        return ToxKind.BODY_WEIGHT_CHANGE
    if matrix is not None and matrix.lower() in {"serum", "plasma", "blood"}:
        return ToxKind.SERUM_CHEMISTRY
    return ToxKind.IN_LIFE_OBSERVATION


def bridge_study_effect_observations(
    study: PharmacologyStudy,
    *,
    endpoint_overrides: dict[str, dict[str, Any]] | None = None,
) -> list[EffectObservation]:
    """Bridge study observations of type 'effect' into EffectObservation records.

    Each effect observation resolves its biosample to locate the linked dose event,
    collection timing, and biological matrix.
    """
    from pharmacoweave.pd_models.records import (
        EffectDirection,
        EffectObservation,
        EffectTimeScale,
    )

    overrides = endpoint_overrides or {}
    sample_by_id: dict[str, Biosample] = {s.sample_id: s for s in study.biosamples}
    dose_by_id: dict[str, DoseEvent] = {d.dose_event_id: d for d in study.dose_events}
    assay_by_id: dict[str, Assay] = {a.assay_id: a for a in study.assays}

    results: list[EffectObservation] = []
    for obs in study.observations:
        if obs.observation_type != ObservationType.EFFECT:
            continue

        biosample = sample_by_id.get(obs.sample_id)
        if biosample is None:
            raise IngestError(
                f"Effect observation {obs.observation_id} references unknown sample_id {obs.sample_id}."
            )
        dose_event = dose_by_id.get(biosample.dose_event_id)
        if dose_event is None:
            raise IngestError(
                f"Biosample {biosample.sample_id} for effect {obs.observation_id} "
                f"references unknown dose_event_id {biosample.dose_event_id}."
            )
        assay = assay_by_id.get(obs.assay_id or "")

        # Look up overrides by observation_id, assay_id, or fallback
        override: dict[str, Any] = (
            overrides.get(obs.observation_id)
            or (overrides.get(obs.assay_id or "") if obs.assay_id else None)
            or (overrides.get(assay.assay_id) if assay else None)
            or {}
        )

        endpoint_name = (
            override.get("endpoint_name")
            or (assay.assay_id if assay else None)
            or obs.raw_value
            or f"effect_{obs.observation_id}"
        )
        matrix = override.get("matrix") or biosample.matrix
        method = (
            override.get("method")
            or (assay.method if assay else None)
            or "in_vivo_bioassay"
        )
        time_scale = override.get("time_scale") or EffectTimeScale.DIRECT
        direction = override.get("direction") or EffectDirection.UNSPECIFIED
        dimension = override.get("dimension") or infer_effect_dimension(obs.unit)
        exposure_ref = override.get("exposure_ref") or biosample.dose_event_id

        results.append(
            EffectObservation(
                effect_id=obs.observation_id,
                endpoint_name=endpoint_name,
                matrix=matrix,
                method=method,
                time_scale=time_scale,
                dimension=dimension,
                direction=direction,
                status=obs.status,
                exposure_ref=exposure_ref,
                limitations=[RESEARCH_USE_LIMITATION],
                numeric_value=obs.numeric_value,
                unit=obs.unit,
                sample_id=obs.sample_id,
                time_since_exposure=biosample.time_since_dose,
                time_since_exposure_unit=biosample.time_since_dose_unit,
                exclusion_reason=obs.exclusion_reason,
            )
        )
    return results


def bridge_study_toxicology_endpoints(
    study: PharmacologyStudy,
    *,
    endpoint_overrides: dict[str, dict[str, Any]] | None = None,
) -> list[ToxEndpoint]:
    """Bridge study observations of type 'toxicity' into ToxEndpoint records.

    Each toxicity observation resolves its biosample to locate the subject and
    compound administered in the linked dose event.
    """
    from pharmacoweave.toxicology.records import ToxEndpoint, ToxKind, ToxSeverity

    overrides = endpoint_overrides or {}
    sample_by_id: dict[str, Biosample] = {s.sample_id: s for s in study.biosamples}
    dose_by_id: dict[str, DoseEvent] = {d.dose_event_id: d for d in study.dose_events}
    assay_by_id: dict[str, Assay] = {a.assay_id: a for a in study.assays}

    results: list[ToxEndpoint] = []
    for obs in study.observations:
        if obs.observation_type != ObservationType.TOXICITY:
            continue

        biosample = sample_by_id.get(obs.sample_id)
        if biosample is None:
            raise IngestError(
                f"Toxicity observation {obs.observation_id} references unknown sample_id {obs.sample_id}."
            )
        dose_event = dose_by_id.get(biosample.dose_event_id)
        if dose_event is None:
            raise IngestError(
                f"Biosample {biosample.sample_id} for toxicity {obs.observation_id} "
                f"references unknown dose_event_id {biosample.dose_event_id}."
            )
        assay = assay_by_id.get(obs.assay_id or "")

        override: dict[str, Any] = (
            overrides.get(obs.observation_id)
            or (overrides.get(obs.assay_id or "") if obs.assay_id else None)
            or (overrides.get(assay.assay_id) if assay else None)
            or {}
        )

        kind = override.get("kind") or infer_tox_kind(biosample.matrix, obs.unit)
        finding = (
            override.get("finding")
            or obs.reported_value
            or (assay.assay_id if assay else None)
            or "toxicology_finding"
        )
        severity = override.get("severity") or (
            ToxSeverity.MINIMAL if kind == ToxKind.HISTOPATHOLOGY else ToxSeverity.NOT_GRADED
        )
        compound_id = override.get("compound_id") or dose_event.compound_id
        subject_id = override.get("subject_id") or biosample.subject_id
        dose_event_id = override.get("dose_event_id") or biosample.dose_event_id

        incidence_affected = override.get("incidence_affected")
        incidence_examined = override.get("incidence_examined")
        if kind == ToxKind.MORTALITY and incidence_affected is None:
            incidence_affected = 0
            incidence_examined = 1

        results.append(
            ToxEndpoint(
                endpoint_id=obs.observation_id,
                compound_id=compound_id,
                kind=kind,
                finding=finding,
                severity=severity,
                limitations=[TOX_RESEARCH_LIMITATION],
                subject_id=subject_id,
                dose_event_id=dose_event_id,
                exposure_ref=override.get("exposure_ref"),
                numeric_value=obs.numeric_value,
                unit=obs.unit,
                incidence_affected=incidence_affected,
                incidence_examined=incidence_examined,
            )
        )
    return results
