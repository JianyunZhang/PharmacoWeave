"""Prior-prediction versus posterior-estimate calibration audit.

A prior is an ``EvidenceRecord`` with ``evidence_type="predicted"`` (TheraWeave
ADMET, CADD, or an exploratory IVIVE). A posterior is an ``estimated`` record
built from NCA or a compartmental fit. The audit reads both and never edits
them: a prediction stays predicted, and an estimate is not replaced by a
prediction.

Per comparison it reports the ratio ``P/O``, the symmetric fold error
``max(P/O, O/P)``, ``log10(P/O)``, and 2-fold / 3-fold flags. Summaries are
grouped by ``measurement_type`` so that an over-predicted CL and an
under-predicted half-life do not cancel in AFE.

Matching is deliberately strict:

- ``clearance`` (CL) never matches ``apparent_clearance`` (CL/F), and
  ``volume`` never matches ``apparent_volume``, because the measurement types
  differ.
- A generic clearance or volume prior does not match a two-compartment
  ``Q``, ``V1``, or ``V2`` estimate unless the prior's ``claim_type`` names
  that parameter.
- A parent-compound prior does not match an estimate for a metabolite analyte.
"""

import math
from collections.abc import Mapping, Sequence
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import AnalyteRole, FiniteFloat, Identifier, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION, EvidenceRecord, EvidenceType
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import BODY_MASS, convert, convert_with_molecular_weight, parse_unit
from pharmacoweave.errors import ContractViolation, UnitError
from pharmacoweave.version import SCHEMA_VERSION

PREDICTION_AUDIT_LIMITATION = (
    "Calibration audit of prior predictions against research-use in vivo estimates. "
    "Fold error describes agreement in this dataset only; it does not validate the prediction method."
)
_COMPARTMENT_CLAIMS = frozenset({"Q", "V1", "V2"})
_UNITY_TOL = 1e-12


class BiasDirection(StrEnum):
    OVER = "over_prediction"
    UNDER = "under_prediction"
    NONE = "none"


class FoldMetrics(Record):
    """Agreement metrics for one predicted/estimated pair."""

    ratio: float = Field(gt=0)
    fold_error: float = Field(ge=1)
    log10_error: FiniteFloat
    within_2_fold: bool
    within_3_fold: bool
    outlier_2_fold: bool
    bias_direction: BiasDirection


class EstimateContext(Record):
    """What a posterior estimate refers to. Built from the study, never guessed from text."""

    evidence_id: str = Field(min_length=1)
    compound_id: str = Field(min_length=1)
    analyte_id: str | None = None
    analyte_role: AnalyteRole | None = None
    subject_id: str | None = None
    body_weight_kg: float | None = Field(default=None, gt=0)


class PredictionAuditEntry(Record):
    """One prior compared against one posterior estimate."""

    prior_evidence_id: str = Field(min_length=1)
    estimate_evidence_id: str = Field(min_length=1)
    prior_source_project: str = Field(min_length=1)
    compound_id: str = Field(min_length=1)
    measurement_type: str = Field(min_length=1)
    prior_claim_type: str = Field(min_length=1)
    estimate_claim_type: str = Field(min_length=1)
    predicted_value_original: float
    predicted_unit_original: str = Field(min_length=1)
    predicted_value: float = Field(gt=0)
    estimated_value: float = Field(gt=0)
    unit: str = Field(min_length=1)
    conversion: Literal["same_unit", "linear", "body_weight", "molecular_weight"]
    ratio: float = Field(gt=0)
    fold_error: float = Field(ge=1)
    log10_error: FiniteFloat
    within_2_fold: bool
    within_3_fold: bool
    outlier_2_fold: bool
    bias_direction: BiasDirection
    subject_id: str | None = None
    analyte_id: str | None = None
    notes: list[str] = Field(default_factory=list)


class SkippedPredictionAudit(Record):
    """A prior (or a prior/estimate pair) that was not compared, with the reason."""

    prior_evidence_id: str = Field(min_length=1)
    reason: str = Field(min_length=1)
    estimate_evidence_id: str | None = None


class PredictionAuditSummary(Record):
    """Agreement summary for one measurement type."""

    measurement_type: str = Field(min_length=1)
    n_comparisons: int = Field(ge=1)
    afe: float = Field(gt=0)
    aafe: float = Field(ge=1)
    fraction_within_2_fold: float = Field(ge=0, le=1)
    fraction_within_3_fold: float = Field(ge=0, le=1)
    n_outlier_2_fold: int = Field(ge=0)


class PredictionAuditReport(Record):
    """Machine-readable audit exported to ``schemas/prediction_audit.schema.json``."""

    audit_id: Identifier
    n_priors: int = Field(ge=0)
    n_estimates: int = Field(ge=0)
    entries: list[PredictionAuditEntry] = Field(default_factory=list)
    skipped: list[SkippedPredictionAudit] = Field(default_factory=list)
    summaries: list[PredictionAuditSummary] = Field(default_factory=list)
    overall_aafe: float | None = Field(default=None, ge=1)
    overall_fraction_within_2_fold: float | None = Field(default=None, ge=0, le=1)
    limitations: list[str] = Field(min_length=1)
    prior_evidence_type: Literal["predicted"] = "predicted"
    estimate_evidence_type: Literal["estimated"] = "estimated"
    research_use_only: Literal[True] = True
    schema_name: Literal["pharmacology-prediction-audit"] = "pharmacology-prediction-audit"
    schema_version: str = SCHEMA_VERSION

    @model_validator(mode="after")
    def versioned_and_consistent(self) -> "PredictionAuditReport":
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {SCHEMA_VERSION}.")
        if RESEARCH_USE_LIMITATION not in self.limitations:
            raise ValueError(f"limitations must include {RESEARCH_USE_LIMITATION!r}.")
        if (self.overall_aafe is None) != (not self.entries):
            raise ValueError("overall_aafe is reported exactly when there is at least one entry.")
        if sum(item.n_comparisons for item in self.summaries) != len(self.entries):
            raise ValueError("summaries must account for every entry.")
        return self


def compute_fold_metrics(predicted: float, estimated: float) -> FoldMetrics:
    """Return ratio, symmetric fold error, log10 error, and 2/3-fold flags for ``predicted`` vs ``estimated``."""

    for label, value in (("predicted", predicted), ("estimated", estimated)):
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be finite and positive for a fold comparison, got {value}.")
    ratio = predicted / estimated
    fold = max(ratio, 1.0 / ratio)
    if math.isclose(ratio, 1.0, rel_tol=_UNITY_TOL, abs_tol=0.0):
        direction = BiasDirection.NONE
    elif ratio > 1.0:
        direction = BiasDirection.OVER
    else:
        direction = BiasDirection.UNDER
    within_2 = fold <= 2.0 + _UNITY_TOL
    return FoldMetrics(
        ratio=ratio,
        fold_error=fold,
        log10_error=math.log10(ratio),
        within_2_fold=within_2,
        within_3_fold=fold <= 3.0 + _UNITY_TOL,
        outlier_2_fold=not within_2,
        bias_direction=direction,
    )


def summarize_audit_entries(entries: Sequence[PredictionAuditEntry]) -> list[PredictionAuditSummary]:
    """Group entries by ``measurement_type`` and compute AFE, AAFE, and 2/3-fold fractions."""

    groups: dict[str, list[PredictionAuditEntry]] = {}
    for entry in entries:
        groups.setdefault(entry.measurement_type, []).append(entry)
    summaries: list[PredictionAuditSummary] = []
    for measurement_type in sorted(groups):
        items = groups[measurement_type]
        n = len(items)
        logs = [item.log10_error for item in items]
        summaries.append(
            PredictionAuditSummary(
                measurement_type=measurement_type,
                n_comparisons=n,
                afe=10.0 ** (sum(logs) / n),
                aafe=10.0 ** (sum(abs(value) for value in logs) / n),
                fraction_within_2_fold=sum(item.within_2_fold for item in items) / n,
                fraction_within_3_fold=sum(item.within_3_fold for item in items) / n,
                n_outlier_2_fold=sum(item.outlier_2_fold for item in items),
            )
        )
    return summaries


def audit_evidence_predictions(
    priors: Sequence[EvidenceRecord],
    estimates: Sequence[EvidenceRecord],
    *,
    audit_id: str = "AUDIT-001",
    contexts: Mapping[str, EstimateContext] | None = None,
    molecular_weights: Mapping[str, tuple[float, str]] | None = None,
) -> PredictionAuditReport:
    """Compare predicted priors with estimated posteriors. Inputs are read, never modified.

    ``contexts`` maps an estimate ``evidence_id`` to its compound, analyte role,
    subject, and body weight. Without a context, an estimate is matched only
    when its ``entity_id`` is itself the compound id. ``molecular_weights``
    maps a compound id to ``(value, unit)`` for explicit mass/molar conversion.
    """

    for record in priors:
        if record.evidence_type is not EvidenceType.PREDICTED:
            raise ContractViolation(
                f"{record.evidence_id}: only evidence_type=predicted records can be audited as priors, "
                f"not {record.evidence_type.value}."
            )
    for record in estimates:
        if record.evidence_type is not EvidenceType.ESTIMATED:
            raise ContractViolation(
                f"{record.evidence_id}: posterior comparators must be evidence_type=estimated, "
                f"not {record.evidence_type.value}."
            )
    context_map = dict(contexts or {})
    mw_map = dict(molecular_weights or {})
    usable_estimates = [
        record for record in estimates if record.numeric_value is not None and record.unit is not None
    ]

    entries: list[PredictionAuditEntry] = []
    skipped: list[SkippedPredictionAudit] = []
    for prior in priors:
        if prior.numeric_value is None or prior.unit is None:
            skipped.append(
                SkippedPredictionAudit(
                    prior_evidence_id=prior.evidence_id,
                    reason="Prior prediction has no numeric_value or unit.",
                )
            )
            continue
        if not math.isfinite(prior.numeric_value) or prior.numeric_value <= 0:
            skipped.append(
                SkippedPredictionAudit(
                    prior_evidence_id=prior.evidence_id,
                    reason="Prior prediction is not a finite positive number.",
                )
            )
            continue
        candidates: list[tuple[EvidenceRecord, EstimateContext]] = []
        for estimate in usable_estimates:
            context = context_map.get(estimate.evidence_id) or EstimateContext(
                evidence_id=estimate.evidence_id, compound_id=estimate.entity_id
            )
            if _matches(prior, estimate, context):
                candidates.append((estimate, context))
        if not candidates:
            skipped.append(
                SkippedPredictionAudit(
                    prior_evidence_id=prior.evidence_id,
                    reason=(
                        f"No estimated {prior.measurement_type} matched compound {prior.entity_id} "
                        f"(claim_type {prior.claim_type}). CL and CL/F, V and V/F, generic and "
                        "compartment-specific parameters, and parent and metabolite are not mixed."
                    ),
                )
            )
            continue
        for estimate, context in candidates:
            assert estimate.numeric_value is not None and estimate.unit is not None
            try:
                converted, conversion = _convert_prior(
                    prior.numeric_value,
                    prior.unit,
                    estimate.unit,
                    body_weight_kg=context.body_weight_kg,
                    molecular_weight=mw_map.get(context.compound_id),
                )
            except UnitError as exc:
                skipped.append(
                    SkippedPredictionAudit(
                        prior_evidence_id=prior.evidence_id,
                        estimate_evidence_id=estimate.evidence_id,
                        reason=f"Units are not comparable: {exc}",
                    )
                )
                continue
            if estimate.numeric_value <= 0 or not math.isfinite(estimate.numeric_value):
                skipped.append(
                    SkippedPredictionAudit(
                        prior_evidence_id=prior.evidence_id,
                        estimate_evidence_id=estimate.evidence_id,
                        reason="Posterior estimate is not a finite positive number.",
                    )
                )
                continue
            metrics = compute_fold_metrics(converted, estimate.numeric_value)
            notes: list[str] = []
            if conversion == "body_weight":
                notes.append(f"Prior was body-weight normalized; scaled by {context.body_weight_kg:g} kg.")
            if conversion == "molecular_weight":
                notes.append("Prior converted between mass and molar units with the compound molecular weight.")
            if prior.claim_type == "ivive_hepatic_plasma_clearance":
                notes.append("IVIVE predicts hepatic clearance only; non-hepatic pathways are not included.")
            entries.append(
                PredictionAuditEntry(
                    prior_evidence_id=prior.evidence_id,
                    estimate_evidence_id=estimate.evidence_id,
                    prior_source_project=prior.source_project.value,
                    compound_id=context.compound_id,
                    measurement_type=prior.measurement_type,
                    prior_claim_type=prior.claim_type,
                    estimate_claim_type=estimate.claim_type,
                    predicted_value_original=prior.numeric_value,
                    predicted_unit_original=prior.unit,
                    predicted_value=converted,
                    estimated_value=estimate.numeric_value,
                    unit=estimate.unit,
                    conversion=conversion,
                    ratio=metrics.ratio,
                    fold_error=metrics.fold_error,
                    log10_error=metrics.log10_error,
                    within_2_fold=metrics.within_2_fold,
                    within_3_fold=metrics.within_3_fold,
                    outlier_2_fold=metrics.outlier_2_fold,
                    bias_direction=metrics.bias_direction,
                    subject_id=context.subject_id or estimate.subject_id,
                    analyte_id=context.analyte_id,
                    notes=notes,
                )
            )

    summaries = summarize_audit_entries(entries)
    overall_aafe: float | None = None
    overall_within_2: float | None = None
    if entries:
        overall_aafe = 10.0 ** (sum(abs(item.log10_error) for item in entries) / len(entries))
        overall_within_2 = sum(item.within_2_fold for item in entries) / len(entries)
    return PredictionAuditReport(
        audit_id=audit_id,
        n_priors=len(priors),
        n_estimates=len(estimates),
        entries=entries,
        skipped=skipped,
        summaries=summaries,
        overall_aafe=overall_aafe,
        overall_fraction_within_2_fold=overall_within_2,
        limitations=[RESEARCH_USE_LIMITATION, PREDICTION_AUDIT_LIMITATION],
    )


def audit_study_predictions(
    study: PharmacologyStudy,
    estimates: Sequence[EvidenceRecord],
    *,
    priors: Sequence[EvidenceRecord] | None = None,
    audit_id: str = "AUDIT-001",
) -> PredictionAuditReport:
    """Audit predicted priors against estimates, resolving compounds from the study.

    ``priors`` defaults to the predicted records in ``study.external_evidence``.
    Compound, analyte role, subject, and body weight are read from the study
    through each estimate's ``derived_from`` links (``analyte:``,
    ``observation:``, ``subject:``). The study is not modified.
    """

    selected = (
        [record for record in study.external_evidence if record.evidence_type is EvidenceType.PREDICTED]
        if priors is None
        else list(priors)
    )
    contexts = build_estimate_contexts(study, estimates)
    molecular_weights = {
        compound.compound_id: (compound.molecular_weight, compound.molecular_weight_unit)
        for compound in study.compounds
        if compound.molecular_weight is not None and compound.molecular_weight_unit is not None
    }
    return audit_evidence_predictions(
        selected,
        estimates,
        audit_id=audit_id,
        contexts=contexts,
        molecular_weights=molecular_weights,
    )


def build_estimate_contexts(
    study: PharmacologyStudy, estimates: Sequence[EvidenceRecord]
) -> dict[str, EstimateContext]:
    """Resolve each estimate's compound, analyte role, subject, and body weight from the study."""

    analytes = {item.analyte_id: item for item in study.analytes}
    observation_analyte = {
        item.observation_id: item.analyte_id for item in study.observations if item.analyte_id is not None
    }
    body_weights: dict[str, float] = {}
    for subject in study.subjects:
        if subject.body_weight is not None:
            body_weights[subject.subject_id] = convert(subject.body_weight.value, subject.body_weight.unit, "kg")
    contexts: dict[str, EstimateContext] = {}
    for record in estimates:
        analyte_id: str | None = None
        subject_id = record.subject_id
        for link in record.derived_from:
            prefix, _, value = link.partition(":")
            if prefix == "analyte" and analyte_id is None:
                analyte_id = value
            elif prefix == "observation" and analyte_id is None:
                analyte_id = observation_analyte.get(value)
            elif prefix == "subject" and subject_id is None:
                subject_id = value
        analyte = analytes.get(analyte_id) if analyte_id is not None else None
        if analyte is None or analyte.compound_id is None:
            continue
        contexts[record.evidence_id] = EstimateContext(
            evidence_id=record.evidence_id,
            compound_id=analyte.compound_id,
            analyte_id=analyte.analyte_id,
            analyte_role=analyte.analyte_role,
            subject_id=subject_id,
            body_weight_kg=None if subject_id is None else body_weights.get(subject_id),
        )
    return contexts


def _matches(prior: EvidenceRecord, estimate: EvidenceRecord, context: EstimateContext) -> bool:
    if prior.entity_id != context.compound_id:
        return False
    if context.analyte_role is AnalyteRole.METABOLITE:
        return False
    if prior.measurement_type != estimate.measurement_type:
        return False
    if prior.claim_type in _COMPARTMENT_CLAIMS:
        return estimate.claim_type == prior.claim_type
    return estimate.claim_type not in _COMPARTMENT_CLAIMS


def _convert_prior(
    value: float,
    source: str,
    target: str,
    *,
    body_weight_kg: float | None,
    molecular_weight: tuple[float, str] | None,
) -> tuple[float, Literal["same_unit", "linear", "body_weight", "molecular_weight"]]:
    if source == target:
        return value, "same_unit"
    source_dims = parse_unit(source).dimensions
    target_dims = parse_unit(target).dimensions
    if source_dims == target_dims:
        return convert(value, source, target), "linear"
    per_body = tuple(left - right for left, right in zip(target_dims, BODY_MASS, strict=True))
    if source_dims == per_body:
        if body_weight_kg is None:
            raise UnitError(
                f"{source} is body-weight normalized and the subject has no body_weight; it is not assumed."
            )
        return convert(value, source, f"{target}/kg") * body_weight_kg, "body_weight"
    if molecular_weight is None:
        raise UnitError(f"Cannot convert {source} to {target} without an explicit compound molecular weight.")
    mw_value, mw_unit = molecular_weight
    converted = convert_with_molecular_weight(
        value, source, target, molecular_weight=mw_value, molecular_weight_unit=mw_unit
    )
    return converted, "molecular_weight"
