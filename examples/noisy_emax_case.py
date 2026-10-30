"""Synthetic noisy direct-effect and sigmoid Emax curves.

Research use only. The direct curve includes a low-concentration noisy observation
that dips slightly below the fixed baseline E0 = 10 percent, plus one BLQ effect
row with no numeric value. Original-scale damped nonlinear least squares recovers
Emax, EC50, and (for the sigmoid curve) the Hill coefficient without replacing
BLQ with zero.
"""

from pharmacoweave.core.entities import ObservationStatus
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.pd_models.records import (
    EffectDimension,
    EffectDirection,
    EffectObservation,
    EffectTimeScale,
)
from pharmacoweave.pkpd.link import ExposureKind, ExposureValue

E0 = 10.0
TRUE_EMAX = 80.0
TRUE_EC50 = 8.0
DIRECT_RECOVERY_RELATIVE_TOLERANCE = 0.12

DIRECT_NOISY_POINTS: tuple[tuple[float, float], ...] = (
    (0.25, 9.450000),
    (1.0, 20.120000),
    (4.0, 35.180000),
    (8.0, 51.940000),
    (16.0, 61.850000),
    (32.0, 75.420000),
    (64.0, 80.260000),
)

SIGMOID_TRUE_EMAX = 80.0
SIGMOID_TRUE_EC50 = 8.0
SIGMOID_TRUE_HILL = 2.0
SIGMOID_RECOVERY_RELATIVE_TOLERANCE = 0.15

SIGMOID_NOISY_POINTS: tuple[tuple[float, float], ...] = (
    (1.0, 10.650000),
    (2.0, 15.420000),
    (4.0, 24.850000),
    (6.0, 40.100000),
    (8.0, 48.900000),
    (12.0, 66.450000),
    (24.0, 80.950000),
    (48.0, 88.500000),
)


def build_noisy_direct_case() -> tuple[list[EffectObservation], list[ExposureValue]]:
    """Return noisy direct-effect observations (including one BLQ) and concentrations."""

    effects: list[EffectObservation] = []
    exposures: list[ExposureValue] = []
    for index, (concentration, observed_effect) in enumerate(DIRECT_NOISY_POINTS, start=1):
        exposure_id = f"NOISY-CONC-{index}"
        exposures.append(
            ExposureValue(
                source_id=exposure_id,
                kind=ExposureKind.CONCENTRATION,
                value=concentration,
                unit="mg/L",
            )
        )
        effects.append(
            EffectObservation(
                effect_id=f"NOISY-EFF-{index}",
                endpoint_name="synthetic noisy response",
                matrix="plasma",
                method="synthetic direct effect",
                time_scale=EffectTimeScale.DIRECT,
                dimension=EffectDimension.PERCENT,
                direction=EffectDirection.INCREASE,
                status=ObservationStatus.QUANTIFIED,
                exposure_ref=exposure_id,
                limitations=[RESEARCH_USE_LIMITATION],
                numeric_value=observed_effect,
                unit="percent",
                time_since_exposure=0.0,
                time_since_exposure_unit="h",
            )
        )
    effects.append(
        EffectObservation(
            effect_id="NOISY-EFF-BLQ",
            endpoint_name="synthetic noisy response",
            matrix="plasma",
            method="synthetic direct effect",
            time_scale=EffectTimeScale.DIRECT,
            dimension=EffectDimension.PERCENT,
            direction=EffectDirection.INCREASE,
            status=ObservationStatus.BLQ,
            exposure_ref="NOISY-CONC-1",
            limitations=[RESEARCH_USE_LIMITATION],
            time_since_exposure=0.0,
            time_since_exposure_unit="h",
        )
    )
    return effects, exposures


def build_noisy_sigmoid_case() -> tuple[list[EffectObservation], list[ExposureValue]]:
    """Return noisy sigmoid Emax observations and concentrations."""

    effects: list[EffectObservation] = []
    exposures: list[ExposureValue] = []
    for index, (concentration, observed_effect) in enumerate(SIGMOID_NOISY_POINTS, start=1):
        exposure_id = f"SIG-CONC-{index}"
        exposures.append(
            ExposureValue(
                source_id=exposure_id,
                kind=ExposureKind.CONCENTRATION,
                value=concentration,
                unit="mg/L",
            )
        )
        effects.append(
            EffectObservation(
                effect_id=f"SIG-EFF-{index}",
                endpoint_name="synthetic sigmoid response",
                matrix="plasma",
                method="synthetic sigmoid effect",
                time_scale=EffectTimeScale.DIRECT,
                dimension=EffectDimension.PERCENT,
                direction=EffectDirection.INCREASE,
                status=ObservationStatus.QUANTIFIED,
                exposure_ref=exposure_id,
                limitations=[RESEARCH_USE_LIMITATION],
                numeric_value=observed_effect,
                unit="percent",
                time_since_exposure=0.0,
                time_since_exposure_unit="h",
            )
        )
    return effects, exposures
