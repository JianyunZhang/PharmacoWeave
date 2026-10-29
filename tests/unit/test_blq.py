"""BLQ policies keep the source observation intact."""

import pytest
from examples.iv_bolus_case import build_study
from pydantic import ValidationError

from pharmacoweave.bioanalysis import BlqPolicy, BlqRule, apply_blq_policy
from pharmacoweave.core.entities import ObservationStatus


def test_default_policy_does_not_substitute_or_mutate() -> None:
    study = build_study()
    before = [item.model_dump() for item in study.observations]
    policy = BlqPolicy(rule=BlqRule.KEEP_AS_CENSORED, rationale="Leave BLQ censored.")
    derived = apply_blq_policy(list(study.observations), list(study.assays), policy)
    assert [item.model_dump() for item in study.observations] == before
    censored = [item for item in derived if item.rule == BlqRule.KEEP_AS_CENSORED.value]
    assert censored
    assert all(item.analysis_value is None for item in censored)


def test_half_lloq_is_explicit_and_audited() -> None:
    study = build_study()
    policy = BlqPolicy(
        rule=BlqRule.LLOQ_OVER_2,
        rationale="Sensitivity analysis using LLOQ/2.",
        acknowledge_bias=True,
    )
    derived = apply_blq_policy(list(study.observations), list(study.assays), policy)
    substituted = [item for item in derived if item.rule == BlqRule.LLOQ_OVER_2.value]
    assert substituted
    assert {item.analysis_value for item in substituted} == {0.005}
    assert all(item.caution for item in substituted)
    stored = [item for item in study.observations if item.status is ObservationStatus.BLQ]
    assert all(item.numeric_value is None for item in stored)


def test_zero_substitution_requires_acknowledgement() -> None:
    with pytest.raises(ValidationError, match="acknowledge_bias"):
        BlqPolicy(rule=BlqRule.ZERO, rationale="Replace BLQ with zero.", acknowledge_bias=False)
    study = build_study()
    policy = BlqPolicy(
        rule=BlqRule.ZERO,
        rationale="Explicit zero substitution for a sensitivity check.",
        acknowledge_bias=True,
    )
    derived = apply_blq_policy(list(study.observations), list(study.assays), policy)
    zeros = [item for item in derived if item.rule == BlqRule.ZERO.value]
    assert zeros and all(item.analysis_value == 0.0 for item in zeros)
    assert all(item.numeric_value is None for item in study.observations if item.status is ObservationStatus.BLQ)
