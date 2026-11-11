"""Weighted two-stage summary. The unweighted summary stays unchanged."""

import json
import math
from pathlib import Path

import pytest

from pharmacoweave.errors import AnalysisError
from pharmacoweave.population_pk.two_stage import (
    IndividualParameterEstimate,
    summarize_log_parameters,
    summarize_weighted_log_parameters,
)

_REFERENCE = json.loads(
    (Path(__file__).resolve().parents[2] / "tests" / "regression" / "weighted_two_stage.json").read_text(
        encoding="utf-8"
    )
)


def test_equal_relative_standard_errors_match_the_unweighted_summary() -> None:
    case = _REFERENCE["equal_relative_standard_error"]
    estimates = [
        _estimate(f"SUBJ-0{index}", value, error)
        for index, (value, error) in enumerate(
            zip(case["estimates_l_h"], case["standard_errors_l_h"], strict=True),
            start=1,
        )
    ]
    weighted = summarize_weighted_log_parameters(estimates)[0]
    unweighted = summarize_log_parameters(estimates)[0]
    hand_variance = 2.0 * (math.log(2.0) ** 2)
    assert weighted.geometric_mean == pytest.approx(case["geometric_mean_l_h"], rel=1e-12, abs=1e-12)
    assert weighted.log_variance == pytest.approx(hand_variance, rel=1e-12, abs=1e-12)
    assert weighted.geometric_mean == pytest.approx(unweighted.geometric_mean, rel=1e-12, abs=1e-12)
    assert weighted.log_variance == pytest.approx(unweighted.log_variance, rel=1e-12, abs=1e-12)
    assert weighted.n_subjects == 2
    assert any("not an NLME posterior" in item for item in weighted.limitations)

    scaled = [
        _estimate(f"SUBJ-0{index}", value, error * 5.0)
        for index, (value, error) in enumerate(
            zip(case["estimates_l_h"], case["standard_errors_l_h"], strict=True),
            start=1,
        )
    ]
    same_weights = summarize_weighted_log_parameters(scaled)[0]
    assert same_weights.geometric_mean == pytest.approx(weighted.geometric_mean, rel=1e-12, abs=1e-12)
    assert same_weights.log_variance == pytest.approx(weighted.log_variance, rel=1e-12, abs=1e-12)


def test_unequal_standard_errors_match_the_hand_formula() -> None:
    pair = _REFERENCE["two_subject_unequal_standard_error"]
    pair_summary = summarize_weighted_log_parameters(
        [
            _estimate(f"SUBJ-0{index}", value, error)
            for index, (value, error) in enumerate(
                zip(pair["estimates_l_h"], pair["standard_errors_l_h"], strict=True),
                start=1,
            )
        ]
    )[0]
    pair_mean = _weighted_geometric_mean(pair["estimates_l_h"], pair["weights"])
    pair_variance = (math.log(pair["estimates_l_h"][1]) - math.log(pair["estimates_l_h"][0])) ** 2 / 2.0
    assert pair_summary.geometric_mean == pytest.approx(pair_mean, rel=1e-12, abs=1e-12)
    assert pair_summary.log_variance == pytest.approx(pair_variance, rel=1e-12, abs=1e-12)
    assert pair_summary.geometric_mean != pytest.approx(0.10, rel=1e-6)

    triple = _REFERENCE["three_subject_unequal_standard_error"]
    triple_summary = summarize_weighted_log_parameters(
        [
            _estimate(f"SUBJ-0{index}", value, error)
            for index, (value, error) in enumerate(
                zip(triple["estimates_l_h"], triple["standard_errors_l_h"], strict=True), start=1
            )
        ]
    )[0]
    assert triple_summary.geometric_mean == pytest.approx(0.1 * math.sqrt(2.0), rel=1e-12, abs=1e-12)
    assert triple_summary.log_variance == pytest.approx((7.0 / 6.0) * (math.log(2.0) ** 2), rel=1e-12, abs=1e-12)
    assert triple_summary.n_subjects == 3
    unweighted = summarize_log_parameters(
        [_estimate(f"SUBJ-0{index}", value, None) for index, value in enumerate(triple["estimates_l_h"], start=1)]
    )[0]
    assert triple_summary.geometric_mean != pytest.approx(unweighted.geometric_mean, rel=1e-6)
    assert triple_summary.log_variance != pytest.approx(unweighted.log_variance, rel=1e-6)


def test_weighted_summary_converts_standard_errors_with_the_estimate() -> None:
    converted = summarize_weighted_log_parameters(
        [
            _estimate("SUBJ-01", 50.0, 10.0, unit="mL/h"),
            _estimate("SUBJ-02", 200.0, 40.0, unit="mL/h"),
        ]
    )[0]
    assert converted.unit == "L/h"
    assert converted.geometric_mean == pytest.approx(0.10, rel=1e-12, abs=1e-12)
    assert converted.log_variance == pytest.approx(2.0 * (math.log(2.0) ** 2), rel=1e-12, abs=1e-12)


def test_missing_or_bad_standard_error_rejects_the_whole_call() -> None:
    valid = [
        _estimate("SUBJ-01", 0.05, 0.01),
        _estimate("SUBJ-02", 0.20, 0.01),
    ]
    with pytest.raises(AnalysisError, match="positive standard error"):
        summarize_weighted_log_parameters(
            [
                *valid,
                _estimate("SUBJ-01", 1.0, 0.1, name="V", unit="L"),
                _estimate("SUBJ-02", 2.0, None, name="V", unit="L"),
            ]
        )
    assert summarize_log_parameters(valid)[0].geometric_mean == pytest.approx(0.10, rel=1e-12, abs=1e-12)
    with pytest.raises(AnalysisError, match="positive standard error"):
        summarize_weighted_log_parameters([valid[0], _estimate("SUBJ-02", 0.20, None)])
    with pytest.raises(AnalysisError, match="positive standard error"):
        summarize_weighted_log_parameters([valid[0], _estimate("SUBJ-02", 0.20, 0.0)])
    with pytest.raises(AnalysisError, match="positive standard error"):
        summarize_weighted_log_parameters([valid[0], _estimate("SUBJ-02", 0.20, -0.01)])
    with pytest.raises(AnalysisError, match="cannot be summarized"):
        summarize_weighted_log_parameters([valid[0], _estimate("SUBJ-02", 0.20, 0.01, name="CL/F")])
    with pytest.raises(AnalysisError, match="fewer than two"):
        summarize_weighted_log_parameters([valid[0]])
    with pytest.raises(AnalysisError, match="repeats"):
        summarize_weighted_log_parameters([valid[0], _estimate("SUBJ-01", 0.20, 0.01)])


def _estimate(
    subject_id: str,
    value: float,
    standard_error: float | None,
    *,
    name: str = "CL",
    unit: str = "L/h",
) -> IndividualParameterEstimate:
    return IndividualParameterEstimate(
        subject_id=subject_id,
        parameter_name=name,  # type: ignore[arg-type]
        estimate=value,
        unit=unit,
        standard_error=standard_error,
    )


def _weighted_geometric_mean(values: list[float], weights: list[float]) -> float:
    logs = [math.log(value) for value in values]
    total = math.fsum(weights)
    mean_log = math.fsum(weight * log for weight, log in zip(weights, logs, strict=True)) / total
    return math.exp(mean_log)
