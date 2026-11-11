"""Closed form for linear turnover with a constant inhibition of production.

    dR/dt = kin * (1 - I) - kout * R

``I = 0`` is the uninhibited linear turnover equation. Parameters are not estimated.
"""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Record
from pharmacoweave.errors import AnalysisError
from pharmacoweave.qsp.turnover import TurnoverPoint


class InhibitedTurnoverSimulation(Record):
    """A predicted turnover curve with constant production inhibition."""

    context_of_use: str = Field(min_length=1)
    evidence_grade: Literal["predicted"] = "predicted"
    kin: float
    kout: float
    r0: float
    inhibition: float
    steady_state: float
    curve: list[TurnoverPoint]
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def simulate_constant_inhibition_turnover(
    *,
    kin: float,
    kout: float,
    r0: float,
    inhibition: float,
    times_h: list[float],
    context_of_use: str,
    estimate_parameters: bool = False,
) -> InhibitedTurnoverSimulation:
    """Evaluate constant-inhibition turnover. Estimation requests are rejected."""

    if estimate_parameters:
        raise AnalysisError("Turnover parameters are inputs. This function does not estimate them.")
    if not math.isfinite(kin) or kin <= 0 or not math.isfinite(kout) or kout <= 0:
        raise AnalysisError("kin and kout must be positive.")
    if not math.isfinite(r0) or r0 < 0:
        raise AnalysisError("R0 must be zero or positive.")
    if not math.isfinite(inhibition) or inhibition < 0 or inhibition >= 1:
        raise AnalysisError("Inhibition must be in [0, 1).")
    if len(times_h) < 1 or any(not math.isfinite(time) or time < 0 for time in times_h):
        raise AnalysisError("Turnover times must be finite and zero or positive.")
    if any(later <= earlier for earlier, later in zip(times_h, times_h[1:], strict=False)):
        raise AnalysisError("Turnover times must strictly increase.")
    production = kin * (1.0 - inhibition)
    steady = production / kout
    curve = [
        TurnoverPoint(time_h=time, response=_response(r0, production, kout, time))
        for time in times_h
    ]
    return InhibitedTurnoverSimulation(
        context_of_use=context_of_use,
        kin=kin,
        kout=kout,
        r0=r0,
        inhibition=inhibition,
        steady_state=steady,
        curve=curve,
        limitations=[
            "The curve solves dR/dt = kin*(1-I) - kout*R in closed form.",
            "Inhibition is a caller-supplied constant. The function does not read a dose.",
            "The evidence grade is predicted.",
            "Parameters were supplied by the caller and were not estimated.",
            "The curve is not a dosing recommendation and is not linked to a clinical endpoint.",
        ],
    )


def _response(r0: float, production: float, kout: float, time_h: float) -> float:
    decay = math.exp(-kout * time_h)
    return r0 * decay + (production / kout) * (1.0 - decay)
