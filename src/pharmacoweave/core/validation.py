"""Relational checks that a single record cannot make.

The report lists errors and warnings. Strict mode rewrites warnings as
errors. Identical subject identifiers across projects are not treated as a
match; only an accepted identity link is eligible, and this function does
not join on it.
"""

import math
from collections import defaultdict
from enum import StrEnum

from pydantic import Field

from pharmacoweave.core.entities import (
    AssayQuantitation,
    DataClass,
    ObservationStatus,
    ObservationType,
    QcStatus,
    Record,
)
from pharmacoweave.core.evidence import ReviewStatus
from pharmacoweave.core.study import PharmacologyStudy
from pharmacoweave.core.units import AMOUNT, DRUG_MASS, convert, convert_with_molecular_weight, parse_unit
from pharmacoweave.errors import ContractViolation, UnitError
from pharmacoweave.integrations.guards import screen_record

KNOWN_MATRICES = frozenset(
    {
        "plasma",
        "serum",
        "whole_blood",
        "blood",
        "csf",
        "urine",
        "bile",
        "liver",
        "brain",
        "kidney",
        "muscle",
        "heart",
        "lung",
        "skin",
        "other",
    }
)
_CLOCK_TOLERANCE_SECONDS = 1.0


class ValidationMode(StrEnum):
    RESEARCH = "research"
    STRICT = "strict"


class ValidationIssue(Record):
    """One validation finding."""

    severity: str
    code: str
    message: str
    entity: str | None = None
    entity_id: str | None = None

    def rewritten_as_error(self) -> "ValidationIssue":
        if self.severity == "error":
            return self
        return self.model_copy(
            update={
                "severity": "error",
                "message": f"{self.message} Promoted to an error in strict mode.",
            }
        )


class ValidationReport(Record):
    """The outcome of validating one study."""

    mode: ValidationMode
    issues: list[ValidationIssue] = Field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    def summary(self) -> str:
        errors = sum(issue.severity == "error" for issue in self.issues)
        warnings = sum(issue.severity == "warning" for issue in self.issues)
        state = "passed" if self.ok else "failed"
        return f"Validation {state} in {self.mode.value} mode: {errors} errors, {warnings} warnings."

    def __str__(self) -> str:
        lines = [self.summary()]
        lines.extend(
            f"- {issue.severity} [{issue.code}] {issue.entity_id or issue.entity or ''}: {issue.message}"
            for issue in self.issues
        )
        return "\n".join(lines)


def validate_study(study: PharmacologyStudy, *, mode: str = "research") -> ValidationReport:
    """Validate foreign keys, clocks, assay limits, identity links, and evidence."""

    selected = ValidationMode(mode)
    issues: list[ValidationIssue] = []
    _unique_ids(study, issues)
    _study_context(study, issues)
    _foreign_keys(study, issues)
    _doses(study, issues)
    _samples(study, issues)
    _urine_collections(study, issues)
    _observations(study, issues)
    _identity(study, issues)
    _evidence(study, issues)
    if selected is ValidationMode.STRICT:
        issues = [issue.rewritten_as_error() for issue in issues]
    return ValidationReport(mode=selected, issues=issues)


def _unique_ids(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    _duplicates("subject", [item.subject_id for item in study.subjects], issues)
    _duplicates("compound", [item.compound_id for item in study.compounds], issues)
    _duplicates("analyte", [item.analyte_id for item in study.analytes], issues)
    _duplicates("assay", [item.assay_id for item in study.assays], issues)
    _duplicates("dose_event", [item.dose_event_id for item in study.dose_events], issues)
    _duplicates("biosample", [item.sample_id for item in study.biosamples], issues)
    _duplicates("observation", [item.observation_id for item in study.observations], issues)
    _duplicates("urine_collection", [item.collection_id for item in study.urine_collections], issues)
    _duplicates("identity_link", [item.link_id for item in study.identity_links], issues)
    _duplicates("evidence", [item.evidence_id for item in study.external_evidence], issues)


def _duplicates(entity: str, identifiers: list[str], issues: list[ValidationIssue]) -> None:
    seen: set[str] = set()
    for identifier in identifiers:
        if identifier in seen:
            _add(
                issues,
                "error",
                "duplicate_id",
                f"Duplicate {entity} id {identifier}.",
                entity,
                identifier,
            )
        seen.add(identifier)


def _study_context(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    if study.metadata.data_class is DataClass.CLINICAL:
        _add(
            issues,
            "warning",
            "clinical_research_use",
            "Clinical records remain research-use inside PharmacoWeave and are not a submission dataset.",
            "study",
            study.metadata.study_id,
        )
    for subject in study.subjects:
        if subject.study_id != study.metadata.study_id:
            _add(
                issues,
                "error",
                "study_id_mismatch",
                f"Subject {subject.subject_id} is tagged with study_id {subject.study_id}.",
                "subject",
                subject.subject_id,
            )
        if subject.species != study.metadata.species:
            _add(
                issues,
                "warning",
                "species_mismatch",
                f"Subject species {subject.species!r} differs from the study species {study.metadata.species!r}.",
                "subject",
                subject.subject_id,
            )


def _foreign_keys(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    subjects = {item.subject_id for item in study.subjects}
    compounds = {item.compound_id for item in study.compounds}
    analytes = {item.analyte_id for item in study.analytes}
    assays = {item.assay_id for item in study.assays}
    doses = {item.dose_event_id for item in study.dose_events}
    samples = {item.sample_id for item in study.biosamples}
    timebases = {item.timebase_id for item in study.metadata.timebases}
    for analyte in study.analytes:
        if analyte.compound_id is not None and analyte.compound_id not in compounds:
            _missing(issues, "analyte", analyte.analyte_id, "compound_id", analyte.compound_id)
    for assay in study.assays:
        if assay.analyte_id not in analytes:
            _missing(issues, "assay", assay.assay_id, "analyte_id", assay.analyte_id)
    for dose in study.dose_events:
        if dose.subject_id not in subjects:
            _missing(issues, "dose_event", dose.dose_event_id, "subject_id", dose.subject_id)
        if dose.compound_id not in compounds:
            _missing(issues, "dose_event", dose.dose_event_id, "compound_id", dose.compound_id)
    for sample in study.biosamples:
        if sample.subject_id not in subjects:
            _missing(issues, "biosample", sample.sample_id, "subject_id", sample.subject_id)
        if sample.dose_event_id not in doses:
            _missing(issues, "biosample", sample.sample_id, "dose_event_id", sample.dose_event_id)
        if sample.timebase_id not in timebases:
            _missing(issues, "biosample", sample.sample_id, "timebase_id", sample.timebase_id)
        _matrix(issues, "biosample", sample.sample_id, sample.matrix)
    for observation in study.observations:
        if observation.sample_id not in samples:
            _missing(issues, "observation", observation.observation_id, "sample_id", observation.sample_id)
        if observation.analyte_id is not None and observation.analyte_id not in analytes:
            _missing(issues, "observation", observation.observation_id, "analyte_id", observation.analyte_id)
        if observation.assay_id is not None and observation.assay_id not in assays:
            _missing(issues, "observation", observation.observation_id, "assay_id", observation.assay_id)
    for analyte in study.analytes:
        _matrix(issues, "analyte", analyte.analyte_id, analyte.matrix)
    for assay in study.assays:
        _matrix(issues, "assay", assay.assay_id, assay.matrix)


def _doses(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    subjects = {item.subject_id: item for item in study.subjects}
    for dose in study.dose_events:
        if dose.body_weight_normalized_value is None:
            continue
        subject = subjects.get(dose.subject_id)
        if subject is None or subject.body_weight is None:
            _add(
                issues,
                "error",
                "dose_normalized_without_body_weight",
                "A body-weight-normalized dose requires the subject's body weight.",
                "dose_event",
                dose.dose_event_id,
            )
            continue
        body_kg = convert(subject.body_weight.value, subject.body_weight.unit, "kg")
        if parse_unit(dose.dose_unit).dimensions == DRUG_MASS:
            absolute = convert(dose.dose_value, dose.dose_unit, "mg")
            stated = convert(
                dose.body_weight_normalized_value,
                dose.body_weight_normalized_unit or "",
                "mg/kg",
            )
            expected = absolute / body_kg
        elif parse_unit(dose.dose_unit).dimensions == AMOUNT:
            absolute = convert(dose.dose_value, dose.dose_unit, "mol")
            stated = convert(
                dose.body_weight_normalized_value,
                dose.body_weight_normalized_unit or "",
                "mol/kg",
            )
            expected = absolute / body_kg
        else:
            continue
        if not math.isclose(expected, stated, rel_tol=1e-3, abs_tol=1e-6):
            _add(
                issues,
                "error",
                "dose_body_weight_mismatch",
                f"Normalized dose {stated} does not match administered dose divided by body weight ({expected}).",
                "dose_event",
                dose.dose_event_id,
            )


def _samples(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    doses = {item.dose_event_id: item for item in study.dose_events}
    for sample in study.biosamples:
        dose = doses.get(sample.dose_event_id)
        if dose is None:
            continue
        if dose.subject_id != sample.subject_id:
            _add(
                issues,
                "error",
                "subject_mismatch",
                (
                    f"Sample {sample.sample_id} belongs to {sample.subject_id}, "
                    f"while dose {dose.dose_event_id} belongs to {dose.subject_id}."
                ),
                "biosample",
                sample.sample_id,
            )
            continue
        delta_seconds = (sample.collection_time_absolute - dose.absolute_time).total_seconds()
        relative_seconds = convert(sample.time_since_dose, sample.time_since_dose_unit, "s")
        if abs(delta_seconds - relative_seconds) > _CLOCK_TOLERANCE_SECONDS:
            _add(
                issues,
                "error",
                "sample_time_inconsistent",
                "time_since_dose does not match collection_time_absolute minus the dose time.",
                "biosample",
                sample.sample_id,
            )


def _urine_collections(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    samples = {item.sample_id: item for item in study.biosamples}
    linked: set[str] = set()
    for collection in study.urine_collections:
        if collection.sample_id in linked:
            _add(
                issues,
                "error",
                "urine_collection_duplicate_sample",
                f"Sample {collection.sample_id} has more than one urine collection.",
                "urine_collection",
                collection.collection_id,
            )
        linked.add(collection.sample_id)
        sample = samples.get(collection.sample_id)
        if sample is None:
            _missing(issues, "urine_collection", collection.collection_id, "sample_id", collection.sample_id)
            continue
        if sample.matrix != "urine":
            _add(
                issues,
                "error",
                "urine_collection_matrix",
                f"Urine collection {collection.collection_id} is tied to matrix {sample.matrix}.",
                "urine_collection",
                collection.collection_id,
            )
            continue
        if sample.is_predose:
            _add(
                issues,
                "error",
                "urine_collection_predose",
                "A urine collection interval cannot end before the dose.",
                "urine_collection",
                collection.collection_id,
            )
            continue
        start_h = convert(collection.interval_start_since_dose, collection.interval_start_unit, "h")
        end_h = convert(sample.time_since_dose, sample.time_since_dose_unit, "h")
        if end_h <= start_h:
            _add(
                issues,
                "error",
                "urine_collection_interval",
                "The urine collection interval must end after it starts.",
                "urine_collection",
                collection.collection_id,
            )


def _observations(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    samples = {item.sample_id: item for item in study.biosamples}
    compounds = {item.compound_id: item for item in study.compounds}
    analytes = {item.analyte_id: item for item in study.analytes}
    assays = {item.assay_id: item for item in study.assays}
    grouped: dict[tuple[str, str | None, ObservationType], list[str]] = defaultdict(list)
    for observation in study.observations:
        grouped[(observation.sample_id, observation.analyte_id, observation.observation_type)].append(
            observation.observation_id
        )
        sample = samples.get(observation.sample_id)
        analyte = analytes.get(observation.analyte_id) if observation.analyte_id else None
        assay = assays.get(observation.assay_id) if observation.assay_id else None
        if analyte is not None and sample is not None and analyte.matrix != sample.matrix:
            _add(
                issues,
                "error",
                "matrix_conflict",
                f"Analyte matrix {analyte.matrix} does not match sample matrix {sample.matrix}.",
                "observation",
                observation.observation_id,
            )
        if assay is None:
            continue
        if analyte is not None and assay.analyte_id != analyte.analyte_id:
            _add(
                issues,
                "error",
                "assay_analyte_mismatch",
                f"Assay {assay.assay_id} measures {assay.analyte_id}, not {analyte.analyte_id}.",
                "observation",
                observation.observation_id,
            )
        if analyte is not None and assay.matrix != analyte.matrix:
            _add(
                issues,
                "error",
                "matrix_conflict",
                f"Assay matrix {assay.matrix} does not match analyte matrix {analyte.matrix}.",
                "assay",
                assay.assay_id,
            )
        if assay.qc_status is QcStatus.REJECTED:
            _add(
                issues,
                "warning",
                "assay_rejected",
                f"Observation uses rejected assay {assay.assay_id}.",
                "observation",
                observation.observation_id,
            )
        if (
            observation.observation_type is ObservationType.CONCENTRATION
            and assay.quantitation is not AssayQuantitation.CONCENTRATION
        ):
            _add(
                issues,
                "error",
                "assay_quantitation",
                "A concentration observation requires a concentration assay.",
                "observation",
                observation.observation_id,
            )
        if (
            observation.status is ObservationStatus.QUANTIFIED
            and observation.observation_type is ObservationType.CONCENTRATION
            and observation.numeric_value is not None
            and observation.unit is not None
            and assay.lloq is not None
            and assay.limit_unit is not None
        ):
            compound = compounds.get(analyte.compound_id) if analyte and analyte.compound_id else None
            try:
                measured = convert_with_molecular_weight(
                    observation.numeric_value,
                    observation.unit,
                    assay.limit_unit,
                    molecular_weight=None if compound is None else compound.molecular_weight,
                    molecular_weight_unit=None if compound is None else compound.molecular_weight_unit,
                )
            except UnitError as exc:
                _add(
                    issues,
                    "error",
                    "unit_conversion_failed",
                    str(exc),
                    "observation",
                    observation.observation_id,
                )
                continue
            if measured < assay.lloq and not math.isclose(measured, assay.lloq, rel_tol=1e-6, abs_tol=1e-9):
                _add(
                    issues,
                    "error",
                    "below_lloq",
                    "Quantified value is below the assay LLOQ. Mark the row BLQ instead of storing the number.",
                    "observation",
                    observation.observation_id,
                )
            if (
                assay.uloq is not None
                and measured > assay.uloq
                and not math.isclose(measured, assay.uloq, rel_tol=1e-6, abs_tol=1e-9)
            ):
                _add(
                    issues,
                    "error",
                    "above_uloq",
                    "Quantified value is above the assay ULOQ. Mark the row above_ULOQ.",
                    "observation",
                    observation.observation_id,
                )
    for (sample_id, analyte_id, observation_type), identifiers in grouped.items():
        if len(identifiers) < 2:
            continue
        rows = [item for item in study.observations if item.observation_id in identifiers]
        if any(not item.replicate_id for item in rows):
            _add(
                issues,
                "warning",
                "replicate_unlabeled",
                (
                    f"Repeated {observation_type.value} rows for sample {sample_id} "
                    f"and analyte {analyte_id} need replicate_id."
                ),
                "observation",
                identifiers[0],
            )


def _identity(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    subjects = {item.subject_id for item in study.subjects}
    accepted: dict[tuple[str, str, str], str] = {}
    for link in study.identity_links:
        if link.review_status is ReviewStatus.UNREVIEWED:
            _add(
                issues,
                "warning",
                "identity_unreviewed",
                "An unreviewed identity link is not a match. The same subject text in two projects is not a match.",
                "identity_link",
                link.link_id,
            )
        if link.review_status is not ReviewStatus.ACCEPTED:
            continue
        if link.canonical_subject_id not in subjects:
            _add(
                issues,
                "error",
                "identity_target_missing",
                f"Accepted link points at unknown subject {link.canonical_subject_id}.",
                "identity_link",
                link.link_id,
            )
        key = (link.source_project.value, link.source_study_id, link.source_subject_id)
        previous = accepted.get(key)
        if previous is not None and previous != link.canonical_subject_id:
            _add(
                issues,
                "error",
                "identity_conflict",
                (
                    f"Source subject {link.source_subject_id} is accepted as both "
                    f"{previous} and {link.canonical_subject_id}."
                ),
                "identity_link",
                link.link_id,
            )
        accepted[key] = link.canonical_subject_id


def _evidence(study: PharmacologyStudy, issues: list[ValidationIssue]) -> None:
    for record in study.external_evidence:
        try:
            screen_record(record)
        except ContractViolation as exc:
            _add(issues, "error", "contract_violation", str(exc), "evidence", record.evidence_id)


def _matrix(issues: list[ValidationIssue], entity: str, entity_id: str, matrix: str) -> None:
    if matrix not in KNOWN_MATRICES:
        _add(
            issues,
            "warning",
            "matrix_unknown",
            f"Matrix {matrix!r} is outside the recommended vocabulary.",
            entity,
            entity_id,
        )


def _missing(
    issues: list[ValidationIssue],
    entity: str,
    entity_id: str,
    field: str,
    target: str,
) -> None:
    _add(
        issues,
        "error",
        "missing_foreign_key",
        f"{field} {target} does not exist.",
        entity,
        entity_id,
    )


def _add(
    issues: list[ValidationIssue],
    severity: str,
    code: str,
    message: str,
    entity: str,
    entity_id: str,
) -> None:
    issues.append(
        ValidationIssue(
            severity=severity,
            code=code,
            message=message,
            entity=entity,
            entity_id=entity_id,
        )
    )
