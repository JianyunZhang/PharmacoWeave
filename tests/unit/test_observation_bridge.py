"""Unit tests for bridging study observations to EffectObservation and ToxEndpoint."""

import pytest
from examples.iv_bolus_case import build_study

from pharmacoweave.core.entities import (
    Observation,
    ObservationStatus,
    ObservationType,
)
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION, TOX_RESEARCH_LIMITATION
from pharmacoweave.core.observation_bridge import (
    bridge_study_effect_observations,
    infer_effect_dimension,
    infer_tox_kind,
)
from pharmacoweave.errors import IngestError, UnitError
from pharmacoweave.pd_models.records import EffectDimension, EffectDirection, EffectTimeScale
from pharmacoweave.toxicology.records import ToxKind, ToxSeverity


def test_infer_effect_dimension() -> None:
    assert infer_effect_dimension(None) == EffectDimension.FRACTION
    assert infer_effect_dimension("%") == EffectDimension.PERCENT
    assert infer_effect_dimension("percent") == EffectDimension.PERCENT
    assert infer_effect_dimension("fraction") == EffectDimension.FRACTION
    assert infer_effect_dimension("ratio") == EffectDimension.FRACTION
    assert infer_effect_dimension("ng/mL") == EffectDimension.CONCENTRATION
    assert infer_effect_dimension("umol/L") == EffectDimension.CONCENTRATION
    assert infer_effect_dimension("1/h") == EffectDimension.RATE
    assert infer_effect_dimension("h") == EffectDimension.TIME

    with pytest.raises(UnitError):
        infer_effect_dimension("invalid_unit_xyz")


def test_infer_tox_kind() -> None:
    assert infer_tox_kind("plasma", "mg/dL") == ToxKind.SERUM_CHEMISTRY
    assert infer_tox_kind("serum", "umol/L") == ToxKind.SERUM_CHEMISTRY
    assert infer_tox_kind("whole_body", "kg") == ToxKind.BODY_WEIGHT_CHANGE
    assert infer_tox_kind("whole_body", "percent") == ToxKind.BODY_WEIGHT_CHANGE
    assert infer_tox_kind("cage", "count") == ToxKind.IN_LIFE_OBSERVATION


def test_bridge_effect_observations() -> None:
    study = build_study()
    # Add an effect observation to study.observations using sample SUBJ-01-T0005
    effect_obs = Observation(
        observation_id="EFF_001",
        sample_id="SUBJ-01-T0005",
        observation_type=ObservationType.EFFECT,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=45.5,
        unit="%",
        assay_id="ASSAY-LCMS-01",
    )
    study_with_effect = study.model_copy(
        update={"observations": [*study.observations, effect_obs]}
    )

    effects = bridge_study_effect_observations(study_with_effect)
    assert len(effects) == 1
    eff = effects[0]
    assert eff.effect_id == "EFF_001"
    assert eff.numeric_value == 45.5
    assert eff.unit == "%"
    assert eff.dimension == EffectDimension.PERCENT
    assert eff.time_scale == EffectTimeScale.DIRECT
    assert eff.exposure_ref == "DOSE-SUBJ-01"
    assert eff.time_since_exposure == 5.0
    assert eff.time_since_exposure_unit == "min"
    assert RESEARCH_USE_LIMITATION in eff.limitations
    assert eff.evidence_grade == "observed"


def test_bridge_effect_observations_with_overrides() -> None:
    study = build_study()
    effect_obs = Observation(
        observation_id="EFF_002",
        sample_id="SUBJ-01-T0005",
        observation_type=ObservationType.EFFECT,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=0.75,
        unit="ratio",
    )
    study_with_effect = study.model_copy(
        update={"observations": [*study.observations, effect_obs]}
    )

    effects = study_with_effect.bridge_effects(
        endpoint_overrides={
            "EFF_002": {
                "endpoint_name": "custom_receptor_occupancy",
                "direction": EffectDirection.INCREASE,
            }
        }
    )
    assert len(effects) == 1
    assert effects[0].endpoint_name == "custom_receptor_occupancy"
    assert effects[0].direction == EffectDirection.INCREASE
    assert effects[0].dimension == EffectDimension.FRACTION


def test_bridge_toxicology_endpoints() -> None:
    study = build_study()
    tox_obs = Observation(
        observation_id="TOX_001",
        sample_id="SUBJ-01-T0015",
        observation_type=ObservationType.TOXICITY,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=85.0,
        unit="mg/dL",
        reported_value="Elevated biomarker",
    )
    study_with_tox = study.model_copy(
        update={"observations": [*study.observations, tox_obs]}
    )

    tox_list = study_with_tox.bridge_toxicology(
        endpoint_overrides={
            "TOX_001": {
                "kind": ToxKind.SERUM_CHEMISTRY,
                "finding": "Elevated biomarker",
                "severity": ToxSeverity.MILD,
            }
        }
    )
    assert len(tox_list) == 1
    t = tox_list[0]
    assert t.endpoint_id == "TOX_001"
    assert t.kind == ToxKind.SERUM_CHEMISTRY
    assert t.finding == "Elevated biomarker"
    assert t.severity == ToxSeverity.MILD
    assert t.numeric_value == 85.0
    assert t.unit == "mg/dL"
    assert t.dose_event_id == "DOSE-SUBJ-01"
    assert t.compound_id == "PW-0001"
    assert TOX_RESEARCH_LIMITATION in t.limitations
    assert t.evidence_grade == "observed"


def test_bridge_effect_missing_sample_fails() -> None:
    study = build_study()
    bad_obs = Observation(
        observation_id="EFF_BAD",
        sample_id="NON_EXISTENT_SAMPLE",
        observation_type=ObservationType.EFFECT,
        status=ObservationStatus.QUANTIFIED,
        numeric_value=12.0,
        unit="percent",
    )
    study_bad = study.model_copy(update={"observations": [bad_obs]})
    with pytest.raises(IngestError, match="unknown sample_id"):
        bridge_study_effect_observations(study_bad)
