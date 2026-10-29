"""Auditable BLQ handling.

Applying a policy returns a new analysis column. It does not rewrite
``numeric_value`` or ``status`` on the observation. The default NCA rule,
``keep_as_censored``, leaves the analysis value empty.
"""

from dataclasses import dataclass
from enum import StrEnum

from pydantic import Field, model_validator

from pharmacoweave.core.entities import Assay, Observation, ObservationStatus, Record
from pharmacoweave.errors import PharmacoWeaveError


class BlqRule(StrEnum):
    EXCLUDE = "exclude"
    KEEP_AS_CENSORED = "keep_as_censored"
    LLOQ_OVER_2 = "lloq_over_2"
    LLOQ = "lloq"
    ZERO = "zero"


_BIASED_RULES = frozenset({BlqRule.LLOQ_OVER_2, BlqRule.LLOQ, BlqRule.ZERO})


class BlqPolicy(Record):
    """A named substitution rule. Biased numeric rules must be acknowledged."""

    rule: BlqRule
    rationale: str = Field(min_length=8)
    acknowledge_bias: bool = False

    @model_validator(mode="after")
    def biased_rules_are_explicit(self) -> "BlqPolicy":
        if self.rule in _BIASED_RULES and not self.acknowledge_bias:
            raise ValueError("Numeric BLQ substitution requires acknowledge_bias=true.")
        return self


@dataclass(frozen=True)
class AnalysisObservation:
    """A derived analysis value. The source observation stays unchanged."""

    observation_id: str
    included: bool
    analysis_value: float | None
    analysis_unit: str | None
    rule: str
    caution: str | None


def apply_blq_policy(
    observations: list[Observation],
    assays: list[Assay],
    policy: BlqPolicy,
) -> list[AnalysisObservation]:
    """Derive analysis values without editing the observations."""

    by_assay = {assay.assay_id: assay for assay in assays}
    return [_apply_one(observation, by_assay.get(observation.assay_id or ""), policy) for observation in observations]


def _apply_one(observation: Observation, assay: Assay | None, policy: BlqPolicy) -> AnalysisObservation:
    if observation.status is ObservationStatus.QUANTIFIED:
        return AnalysisObservation(
            observation_id=observation.observation_id,
            included=True,
            analysis_value=observation.numeric_value,
            analysis_unit=observation.unit,
            rule="observed",
            caution=None,
        )
    if observation.status is ObservationStatus.ABOVE_ULOQ:
        return AnalysisObservation(
            observation_id=observation.observation_id,
            included=True,
            analysis_value=None,
            analysis_unit=observation.unit,
            rule="above_uloq",
            caution="Value is above ULOQ. This policy does not substitute a number.",
        )
    if observation.status in {ObservationStatus.MISSING, ObservationStatus.INVALID}:
        return AnalysisObservation(
            observation_id=observation.observation_id,
            included=False,
            analysis_value=None,
            analysis_unit=None,
            rule=observation.status.value,
            caution=observation.exclusion_reason,
        )
    if policy.rule is BlqRule.EXCLUDE:
        return _censored(observation, included=False, rule=policy.rule.value, value=None, unit=None, caution=None)
    if policy.rule is BlqRule.KEEP_AS_CENSORED:
        return _censored(
            observation,
            included=True,
            rule=policy.rule.value,
            value=None,
            unit=observation.unit,
            caution="Censored below LLOQ. No number was substituted.",
        )
    if assay is None or assay.lloq is None or assay.limit_unit is None:
        raise PharmacoWeaveError(f"{observation.observation_id} cannot use {policy.rule.value} without an assay LLOQ.")
    substituted = assay.lloq / 2 if policy.rule is BlqRule.LLOQ_OVER_2 else assay.lloq
    if policy.rule is BlqRule.ZERO:
        substituted = 0.0
    return _censored(
        observation,
        included=True,
        rule=policy.rule.value,
        value=substituted,
        unit=assay.limit_unit,
        caution="Biased numeric substitution for a censored observation. The stored status is still BLQ.",
    )


def _censored(
    observation: Observation,
    *,
    included: bool,
    rule: str,
    value: float | None,
    unit: str | None,
    caution: str | None,
) -> AnalysisObservation:
    return AnalysisObservation(
        observation_id=observation.observation_id,
        included=included,
        analysis_value=value,
        analysis_unit=unit,
        rule=rule,
        caution=caution,
    )
