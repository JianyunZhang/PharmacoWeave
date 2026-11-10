"""Closed form for one linear turnover equation.

    dR/dt = kin - kout * R
    R(t) = R0 * exp(-kout * t) + (kin / kout) * (1 - exp(-kout * t))

The curve is predicted. Parameters are not estimated.
"""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Record
from pharmacoweave.errors import AnalysisError


class TurnoverPoint(Record):
    """One predicted response."""

    time_h: float
    response: float


class TurnoverSimulation(Record):
    """A predicted turnover curve and its steady state kin/kout."""

    context_of_use: str = Field(min_length=1)
    evidence_grade: Literal["predicted"] = "predicted"
    kin: float
    kout: float
    r0: float
    steady_state: float
    curve: list[TurnoverPoint]
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def simulate_linear_turnover(
    *,
    kin: float,
    kout: float,
    r0: float,
    times_h: list[float],
    context_of_use: str,
    estimate_parameters: bool = False,
) -> TurnoverSimulation:
    """Evaluate the linear turnover closed form. Estimation requests are rejected."""

    if estimate_parameters:
        raise AnalysisError("Turnover parameters are inputs. This function does not estimate them.")
    if not math.isfinite(kin) or kin <= 0 or not math.isfinite(kout) or kout <= 0:
        raise AnalysisError("kin and kout must be positive.")
    if not math.isfinite(r0) or r0 < 0:
        raise AnalysisError("R0 must be zero or positive.")
    if len(times_h) < 1 or any(not math.isfinite(time) or time < 0 for time in times_h):
        raise AnalysisError("Turnover times must be finite and zero or positive.")
    if any(later <= earlier for earlier, later in zip(times_h, times_h[1:], strict=False)):
        raise AnalysisError("Turnover times must strictly increase.")
    steady = kin / kout
    curve = [
        TurnoverPoint(time_h=time, response=_response(r0, kin, kout, time))
        for time in times_h
    ]
    return TurnoverSimulation(
        context_of_use=context_of_use,
        kin=kin,
        kout=kout,
        r0=r0,
        steady_state=steady,
        curve=curve,
        limitations=[
            "The curve solves dR/dt = kin - kout*R in closed form.",
            "The evidence grade is predicted.",
            "Parameters were supplied by the caller and were not estimated.",
            "The curve is not a dosing recommendation and is not linked to a clinical endpoint.",
        ],
    )


def _response(r0: float, kin: float, kout: float, time_h: float) -> float:
    decay = math.exp(-kout * time_h)
    return r0 * decay + (kin / kout) * (1.0 - decay)
