"""Effect time scale, dimension, and the observed-versus-predicted boundary."""

import pytest
from pydantic import ValidationError

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.errors import AnalysisError, ContractViolation
from pharmacoweave.pd_models import (
    EffectDimension,
    EffectDirection,
    EffectObservation,
    EffectTimeScale,
    PredictedEffect,
    effect_observation_from_prediction,
    estimate_hill,
)


def _effect(**overrides: object) -> EffectObservation:
    payload: dict[str, object] = {
        "effect_id": "EFF-1",
        "endpoint_name": "response",
        "matrix": "plasma",
        "method": "synthetic",
        "time_scale": EffectTimeScale.DIRECT,
        "dimension": EffectDimension.PERCENT,
        "direction": EffectDirection.INCREASE,
        "status": ObservationStatus.QUANTIFIED,
        "exposure_ref": "CONC-1",
        "limitations": [RESEARCH_USE_LIMITATION],
        "numeric_value": 50.0,
        "unit": "percent",
        "time_since_exposure": 0.0,
        "time_since_exposure_unit": "h",
    }
    payload.update(overrides)
    return EffectObservation(**payload)  # type: ignore[arg-type]


def test_direct_effect_requires_a_time_and_a_matching_unit() -> None:
    effect = _effect()
    assert effect.evidence_grade == "observed"
    assert effect.time_scale is EffectTimeScale.DIRECT
    with pytest.raises(ValidationError, match="time_since_exposure"):
        _effect(time_since_exposure=None, time_since_exposure_unit=None)
    with pytest.raises(ValidationError, match="does not match"):
        _effect(unit="mg/L")


def test_interval_and_time_to_event_keep_their_clocks() -> None:
    interval = _effect(
        time_scale=EffectTimeScale.INTERVAL,
        time_since_exposure=None,
        time_since_exposure_unit=None,
        window_start=0.0,
        window_end=2.0,
        window_unit="h",
    )
    assert interval.window_end == 2.0
    with pytest.raises(ValidationError, match="window_end"):
        _effect(
            time_scale=EffectTimeScale.INTERVAL,
            time_since_exposure=None,
            time_since_exposure_unit=None,
            window_start=2.0,
            window_end=2.0,
            window_unit="h",
        )
    event = _effect(
        time_scale=EffectTimeScale.TIME_TO_EVENT,
        dimension=EffectDimension.TIME,
        numeric_value=3.0,
        unit="h",
        time_since_exposure=None,
        time_since_exposure_unit=None,
    )
    assert event.dimension is EffectDimension.TIME
    with pytest.raises(ValidationError, match="dimension time"):
        _effect(dimension=EffectDimension.TIME)


def test_blq_effect_cannot_store_zero() -> None:
    with pytest.raises(ValidationError, match="numeric_value"):
        _effect(status=ObservationStatus.BLQ, numeric_value=0.0, unit="percent")


def test_prediction_cannot_be_labeled_observed() -> None:
    predicted = PredictedEffect(
        effect_id="EFF-1",
        exposure_ref="CONC-1",
        value=40.0,
        unit="percent",
        derived_from=("fit:direct-emax",),
    )
    assert predicted.evidence_grade == "predicted"
    with pytest.raises(ValueError, match="predicted"):
        PredictedEffect(
            effect_id="EFF-1",
            exposure_ref="CONC-1",
            value=40.0,
            unit="percent",
            derived_from=("fit:direct-emax",),
            evidence_grade="observed",
        )
    with pytest.raises(ContractViolation, match="observed effect"):
        effect_observation_from_prediction(predicted)
    with pytest.raises(ValidationError):
        _effect(evidence_grade="predicted")


def test_hill_estimation_is_refused() -> None:
    with pytest.raises(AnalysisError, match="hill"):
        estimate_hill()
