"""Observed effect records.

An effect is not a concentration and not a model prediction. Time scale and
dimension are required so a percent change, a biomarker concentration, and a
time-to-event stay distinct. The study observation table is unchanged;
these records sit on top of an exposure id from Phase 1.
"""

from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from pharmacoweave.core.entities import FiniteFloat, Identifier, ObservationStatus, Record
from pharmacoweave.core.evidence import RESEARCH_USE_LIMITATION
from pharmacoweave.core.units import CONC_AMOUNT, CONC_MASS, DIMENSIONLESS, RATE, TIME, parse_unit


class EffectTimeScale(StrEnum):
    """How the effect clock relates to the linked exposure."""

    DIRECT = "direct"
    INTERVAL = "interval"
    TIME_TO_EVENT = "time_to_event"


class EffectDimension(StrEnum):
    """The dimension of the effect number. This is not an exposure dimension."""

    FRACTION = "fraction"
    PERCENT = "percent"
    CONCENTRATION = "concentration"
    RATE = "rate"
    TIME = "time"


class EffectDirection(StrEnum):
    INCREASE = "increase"
    DECREASE = "decrease"
    UNSPECIFIED = "unspecified"


def unit_matches_effect_dimension(unit: str, dimension: EffectDimension) -> bool:
    """Return whether ``unit`` is the dimension declared for an effect."""

    parsed = parse_unit(unit)
    if dimension is EffectDimension.FRACTION:
        return parsed.dimensions == DIMENSIONLESS and parsed.factor == 1.0
    if dimension is EffectDimension.PERCENT:
        return parsed.dimensions == DIMENSIONLESS and parsed.factor == 0.01
    if dimension is EffectDimension.CONCENTRATION:
        return parsed.dimensions in {CONC_MASS, CONC_AMOUNT}
    if dimension is EffectDimension.RATE:
        return parsed.dimensions == RATE
    if dimension is EffectDimension.TIME:
        return parsed.dimensions == TIME
    return False


class EffectObservation(Record):
    """One observed effect tied to an exposure id.

    ``evidence_grade`` is fixed to observed. A model prediction is a separate
    object and cannot be stored here.
    """

    effect_id: Identifier
    endpoint_name: str = Field(min_length=1)
    matrix: str = Field(min_length=1)
    method: str = Field(min_length=1)
    time_scale: EffectTimeScale
    dimension: EffectDimension
    direction: EffectDirection
    status: ObservationStatus
    exposure_ref: Identifier
    limitations: list[str] = Field(min_length=1)
    numeric_value: FiniteFloat | None = None
    unit: str | None = None
    sample_id: Identifier | None = None
    time_since_exposure: FiniteFloat | None = None
    time_since_exposure_unit: str | None = None
    window_start: FiniteFloat | None = None
    window_end: FiniteFloat | None = None
    window_unit: str | None = None
    exclusion_reason: str | None = None
    evidence_grade: Literal["observed"] = "observed"
    research_use_only: Literal[True] = True

    @model_validator(mode="after")
    def time_scale_and_value(self) -> "EffectObservation":
        if RESEARCH_USE_LIMITATION not in self.limitations:
            raise ValueError(f"limitations must include {RESEARCH_USE_LIMITATION!r}.")
        _time_scale(self)
        if self.time_scale is not EffectTimeScale.TIME_TO_EVENT and self.dimension is EffectDimension.TIME:
            raise ValueError("dimension time is only valid for time_scale time_to_event.")
        if self.status is ObservationStatus.QUANTIFIED:
            if self.numeric_value is None or self.unit is None:
                raise ValueError("A quantified effect requires numeric_value and unit.")
            if not unit_matches_effect_dimension(self.unit, self.dimension):
                raise ValueError(f"Unit {self.unit} does not match effect dimension {self.dimension.value}.")
            if self.dimension is EffectDimension.CONCENTRATION and self.numeric_value <= 0:
                raise ValueError("A quantified concentration effect must be positive.")
            if self.dimension is EffectDimension.TIME and self.numeric_value <= 0:
                raise ValueError("A quantified time-to-event must be positive.")
        elif self.numeric_value is not None:
            raise ValueError(
                "numeric_value is recorded only when status is quantified. "
                "BLQ, above_ULOQ, missing, and invalid stay non-numeric."
            )
        if self.status is ObservationStatus.INVALID and not self.exclusion_reason:
            raise ValueError("An invalid effect requires exclusion_reason.")
        return self


def _time_scale(record: EffectObservation) -> None:
    window = (record.window_start, record.window_end, record.window_unit)
    clock = (record.time_since_exposure, record.time_since_exposure_unit)
    if record.time_scale is EffectTimeScale.DIRECT:
        if clock[0] is None or clock[1] is None:
            raise ValueError("A direct effect requires time_since_exposure and time_since_exposure_unit.")
        if parse_unit(clock[1]).dimensions != TIME:
            raise ValueError("time_since_exposure_unit must be a time unit.")
        if any(item is not None for item in window):
            raise ValueError("A direct effect does not carry an interval window.")
        return
    if record.time_scale is EffectTimeScale.INTERVAL:
        if any(item is None for item in window):
            raise ValueError("An interval effect requires window_start, window_end, and window_unit.")
        if parse_unit(record.window_unit or "").dimensions != TIME:
            raise ValueError("window_unit must be a time unit.")
        if record.window_end is None or record.window_start is None or record.window_end <= record.window_start:
            raise ValueError("window_end must be greater than window_start.")
        if any(item is not None for item in clock):
            raise ValueError("An interval effect records a window, not time_since_exposure.")
        return
    if record.dimension is not EffectDimension.TIME:
        raise ValueError("time_to_event requires dimension time.")
    if any(item is not None for item in window) or any(item is not None for item in clock):
        raise ValueError("time_to_event stores the duration as the effect value, not as a second clock.")
