"""Auditable replicate collapsing and derived observation table tests."""

import math

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.bioanalysis import BlqPolicy, BlqRule, derive_study_observations
from pharmacoweave.core.entities import AnalyteForm, ObservationStatus
from pharmacoweave.core.profiles import collapse_replicates, iv_bolus_profiles
from pharmacoweave.errors import AnalysisError


def _study_with_replicates():
    study = build_study()
    quantified = [
        obs
        for obs in study.observations
        if obs.status is ObservationStatus.QUANTIFIED and obs.numeric_value is not None
    ]
    target = quantified[0]
    assert target.numeric_value is not None
    rep_b = target.model_copy(
        update={
            "observation_id": f"{target.observation_id}.repB",
            "numeric_value": target.numeric_value * 1.21,
            "replicate_id": "2",
        }
    )
    blq_rep = target.model_copy(
        update={
            "observation_id": f"{target.observation_id}.repBLQ",
            "status": ObservationStatus.BLQ,
            "numeric_value": None,
            "unit": None,
            "raw_value": "<LLOQ",
            "reported_value": "BLQ",
            "replicate_id": "3",
        }
    )
    updated = study.model_copy(update={"observations": [*study.observations, rep_b, blq_rep]})
    return updated, target, rep_b


def test_replicate_policy_error_by_default_and_collapse_arithmetic_or_geometric() -> None:
    study, target, rep_b = _study_with_replicates()
    assert target.numeric_value is not None and rep_b.numeric_value is not None

    with pytest.raises(AnalysisError, match="Collapse replicates"):
        iv_bolus_profiles(study)

    collapsed_arith, summaries = collapse_replicates(study, method="arithmetic_mean")
    assert len(summaries) == 1
    summary = summaries[0]
    expected_arith = 0.5 * (target.numeric_value + rep_b.numeric_value)
    assert summary.collapsed_value == pytest.approx(expected_arith)
    assert summary.n_replicates == 2
    assert summary.n_censored_unchanged == 1
    assert summary.cv_percent > 0.0
    assert summary.source_observation_ids == [target.observation_id, rep_b.observation_id]

    blq_rows = [obs for obs in collapsed_arith.observations if obs.status is ObservationStatus.BLQ]
    assert any(obs.observation_id.endswith(".repBLQ") for obs in blq_rows)
    assert all(obs.numeric_value is None for obs in blq_rows)

    profiles_arith = iv_bolus_profiles(study, replicate_policy="arithmetic_mean")
    assert profiles_arith

    collapsed_geom, geom_summaries = collapse_replicates(study, method="geometric_mean")
    expected_geom = math.sqrt(target.numeric_value * rep_b.numeric_value)
    assert geom_summaries[0].collapsed_value == pytest.approx(expected_geom)
    assert collapsed_geom.check().ok


def test_collapse_replicates_refuses_mixed_analyte_forms() -> None:
    study, target, rep_b = _study_with_replicates()
    mixed_rep = rep_b.model_copy(update={"analyte_form": AnalyteForm.FREE})
    mixed_study = study.model_copy(
        update={
            "observations": [
                obs if obs.observation_id != rep_b.observation_id else mixed_rep
                for obs in study.observations
            ]
        }
    )
    with pytest.raises(AnalysisError, match="analyte_form"):
        collapse_replicates(mixed_study)


def test_derive_study_observations_preserves_raw_and_audits_policy() -> None:
    study = build_study()
    before = [obs.model_dump() for obs in study.observations]

    default_rows = derive_study_observations(study)
    assert [obs.model_dump() for obs in study.observations] == before
    assert len(default_rows) == len(study.observations)

    predose_rows = [row for row in default_rows if row.is_predose]
    assert predose_rows
    assert all(not row.included_in_analysis for row in predose_rows)
    assert all(row.policy_rule.startswith("predose:") for row in predose_rows)

    blq_rows = [row for row in default_rows if row.stored_status is ObservationStatus.BLQ and not row.is_predose]
    assert blq_rows
    assert all(row.stored_numeric_value is None and row.analysis_value is None for row in blq_rows)

    biased = BlqPolicy(
        rule=BlqRule.LLOQ_OVER_2,
        rationale="Sensitivity check for derived observations parquet.",
        acknowledge_bias=True,
    )
    substituted_rows = derive_study_observations(study, policy=biased)
    sub_blq = [row for row in substituted_rows if row.stored_status is ObservationStatus.BLQ and not row.is_predose]
    assert sub_blq
    assert all(row.stored_numeric_value is None for row in sub_blq)
    assert all(row.analysis_value == pytest.approx(0.005) for row in sub_blq)
    assert all(row.caution and "stored status is still BLQ" in row.caution for row in sub_blq)
