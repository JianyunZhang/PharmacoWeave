"""Failure cases for direct Emax and the concentration link."""

import pytest
from examples.direct_effect_and_evidence import E0, build_effects, build_exposures

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.modeling import DiagnosticStatus
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pd_models import (
    DirectEffectPair,
    EffectDimension,
    EffectDirection,
    EffectObservation,
    EffectTimeScale,
    emax_effect,
    fit_direct_emax,
)
from pharmacoweave.pkpd import (
    ExposureKind,
    ExposureValue,
    fit_effect_lag,
    fit_indirect_response,
    fit_linked_direct_emax,
)


def _pair(index: int, concentration: float, effect: float) -> DirectEffectPair:
    return DirectEffectPair(
        effect_id=f"EFF-{index}",
        exposure_ref=f"CONC-{index}",
        concentration=concentration,
        concentration_unit="mg/L",
        effect=effect,
        effect_unit="percent",
    )


def test_sigmoid_forward_equation_at_half_maximal_concentration() -> None:
    assert emax_effect(concentration=10.0, e0=0.0, emax=100.0, ec50=10.0, hill=1.0) == pytest.approx(50.0)
    assert emax_effect(concentration=5.0, e0=0.0, emax=100.0, ec50=10.0, hill=2.0) == pytest.approx(20.0)
    assert emax_effect(concentration=0.0, e0=4.0, emax=100.0, ec50=10.0) == 4.0


def test_flat_effect_and_short_series_do_not_report_ec50() -> None:
    flat = fit_direct_emax(
        [_pair(index, float(index), 20.0) for index in range(1, 5)],
        e0=10.0,
        e0_unit="percent",
        direction=EffectDirection.INCREASE,
    )
    assert flat.diagnostic.status is DiagnosticStatus.FAIL
    assert flat.estimates == ()
    assert flat.predicted == ()
    assert "do not vary" in flat.reason

    short = fit_direct_emax([_pair(1, 1.0, 20.0), _pair(2, 5.0, 40.0)], e0=10.0, e0_unit="percent")
    assert short.diagnostic.status is DiagnosticStatus.FAIL
    assert "Fewer than three" in short.reason
    assert all(item.name != "EC50" for item in short.estimates)


def test_baseline_hit_identical_concentrations_and_wrong_direction_fail() -> None:
    baseline_nls = fit_direct_emax(
        [_pair(1, 1.0, 10.0), _pair(2, 5.0, 40.0), _pair(3, 10.0, 50.0)],
        e0=10.0,
        e0_unit="percent",
    )
    assert baseline_nls.diagnostic.status is DiagnosticStatus.FAIL
    assert baseline_nls.estimates == ()
    assert "plateau" in baseline_nls.reason

    baseline_linearized = fit_direct_emax(
        [_pair(1, 1.0, 10.0), _pair(2, 5.0, 40.0), _pair(3, 10.0, 50.0)],
        e0=10.0,
        e0_unit="percent",
        method="linearized",
    )
    assert baseline_linearized.diagnostic.status is DiagnosticStatus.FAIL
    assert "baseline" in baseline_linearized.reason

    same_concentration = fit_direct_emax(
        [_pair(index, 5.0, 20.0 + index) for index in range(1, 4)],
        e0=10.0,
        e0_unit="percent",
    )
    assert "do not vary" in same_concentration.reason

    decreasing = [
        _pair(index, concentration, emax_effect(concentration=concentration, e0=100.0, emax=-80.0, ec50=8.0))
        for index, concentration in enumerate((1.0, 5.0, 10.0, 20.0), start=1)
    ]
    wrong_direction = fit_direct_emax(
        decreasing,
        e0=100.0,
        e0_unit="percent",
        direction=EffectDirection.INCREASE,
    )
    assert wrong_direction.diagnostic.status is DiagnosticStatus.FAIL
    assert wrong_direction.estimates == ()
    assert "increase" in wrong_direction.reason


def test_blq_is_omitted_and_auc_is_not_a_concentration() -> None:
    effects = build_effects()
    effects.append(
        EffectObservation(
            effect_id="EFF-BLQ",
            endpoint_name="synthetic response",
            matrix="plasma",
            method="synthetic direct effect",
            time_scale=EffectTimeScale.DIRECT,
            dimension=EffectDimension.PERCENT,
            direction=EffectDirection.INCREASE,
            status=ObservationStatus.BLQ,
            exposure_ref="CONC-1",
            limitations=[RESEARCH_USE_LIMITATION],
            time_since_exposure=0.0,
            time_since_exposure_unit="h",
        )
    )
    linked, fit = fit_linked_direct_emax(effects, build_exposures(), e0=E0, e0_unit="percent")
    assert [item.effect_id for item in linked.omitted] == ["EFF-BLQ"]
    assert "not replaced with zero" in linked.omitted[0].reason
    assert fit.diagnostic.status is DiagnosticStatus.PASS
    assert all(point.value != 0.0 for point in fit.predicted)

    auc = ExposureValue(source_id="AUC-1", kind=ExposureKind.AUC, value=5.0, unit="mg*h/L")
    swapped = [effect.model_copy(update={"exposure_ref": "AUC-1"}) for effect in build_effects()]
    with pytest.raises(AnalysisError, match="AUC"):
        fit_linked_direct_emax(swapped, [auc], e0=E0, e0_unit="percent")

    cmax = ExposureValue(source_id="CMAX-1", kind=ExposureKind.CMAX, value=1.0, unit="mg/L")
    swapped_cmax = [effect.model_copy(update={"exposure_ref": "CMAX-1"}) for effect in build_effects()]
    with pytest.raises(AnalysisError, match="Cmax"):
        fit_linked_direct_emax(swapped_cmax, [cmax], e0=E0, e0_unit="percent")


def test_interval_effect_is_not_linked_as_a_direct_concentration() -> None:
    effects = build_effects()
    effects.append(
        EffectObservation(
            effect_id="EFF-WINDOW",
            endpoint_name="synthetic response",
            matrix="plasma",
            method="synthetic direct effect",
            time_scale=EffectTimeScale.INTERVAL,
            dimension=EffectDimension.PERCENT,
            direction=EffectDirection.INCREASE,
            status=ObservationStatus.QUANTIFIED,
            exposure_ref="CONC-1",
            limitations=[RESEARCH_USE_LIMITATION],
            numeric_value=40.0,
            unit="percent",
            window_start=0.0,
            window_end=2.0,
            window_unit="h",
        )
    )
    linked, fit = fit_linked_direct_emax(effects, build_exposures(), e0=E0, e0_unit="percent")
    assert linked.omitted[0].effect_id == "EFF-WINDOW"
    assert "interval" in linked.omitted[0].reason
    assert fit.diagnostic.status is DiagnosticStatus.PASS
    assert all(point.effect_id != "EFF-WINDOW" for point in fit.predicted)


def test_lag_and_indirect_response_are_not_fitted() -> None:
    with pytest.raises(AnalysisError, match="lag"):
        fit_effect_lag()
    with pytest.raises(AnalysisError, match="Indirect-response"):
        fit_indirect_response()
