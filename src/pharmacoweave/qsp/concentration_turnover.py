"""Turnover whose production is multiplied by a one-compartment IV bolus concentration.

    dR/dt = kin * C(t) - kout * R
    C(t) = (Dose / V) * exp(-(CL / V) * t)

Loss-rate modulation by concentration is not evaluated: that integral is not elementary.
Parameters are not estimated.
"""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Record
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import iv_bolus_concentration


class ConcentrationTurnoverPoint(Record):
    """Predicted plasma concentration and turnover response at one time."""

    time_h: float
    plasma_concentration_mg_l: float
    response: float


class ConcentrationTurnoverSimulation(Record):
    """Predicted turnover driven by a decaying IV bolus concentration."""

    context_of_use: str = Field(min_length=1)
    evidence_grade: Literal["predicted"] = "predicted"
    kin: float
    kout: float
    r0: float
    dose_mg: float
    volume_l: float
    clearance_l_h: float
    steady_state: float
    curve: list[ConcentrationTurnoverPoint]
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def simulate_concentration_driven_turnover(
    *,
    kin: float,
    kout: float,
    r0: float,
    dose_mg: float,
    volume_l: float,
    clearance_l_h: float,
    times_h: list[float],
    context_of_use: str,
    estimate_parameters: bool = False,
) -> ConcentrationTurnoverSimulation:
    """Evaluate production driven by the accepted IV bolus curve."""

    if estimate_parameters:
        raise AnalysisError("Turnover parameters are inputs. This function does not estimate them.")
    if not math.isfinite(kin) or kin <= 0 or not math.isfinite(kout) or kout <= 0:
        raise AnalysisError("kin and kout must be positive.")
    if not math.isfinite(r0) or r0 < 0:
        raise AnalysisError("R0 must be zero or positive.")
    if not math.isfinite(dose_mg) or dose_mg <= 0 or not math.isfinite(volume_l) or volume_l <= 0:
        raise AnalysisError("Dose and volume must be positive.")
    if not math.isfinite(clearance_l_h) or clearance_l_h <= 0:
        raise AnalysisError("Clearance must be positive.")
    if len(times_h) < 1 or any(not math.isfinite(time) or time < 0 for time in times_h):
        raise AnalysisError("Turnover times must be finite and zero or positive.")
    if any(later <= earlier for earlier, later in zip(times_h, times_h[1:], strict=False)):
        raise AnalysisError("Turnover times must strictly increase.")
    kel = clearance_l_h / volume_l
    amount_over_volume = dose_mg / volume_l
    curve = [
        ConcentrationTurnoverPoint(
            time_h=time,
            plasma_concentration_mg_l=iv_bolus_concentration(
                dose=dose_mg,
                dose_unit="mg",
                volume=volume_l,
                volume_unit="L",
                clearance=clearance_l_h,
                clearance_unit="L/h",
                time=time,
                time_unit="h",
                concentration_unit="mg/L",
            ),
            response=_response(r0, kin, kout, amount_over_volume, kel, time),
        )
        for time in times_h
    ]
    return ConcentrationTurnoverSimulation(
        context_of_use=context_of_use,
        kin=kin,
        kout=kout,
        r0=r0,
        dose_mg=dose_mg,
        volume_l=volume_l,
        clearance_l_h=clearance_l_h,
        steady_state=0.0,
        curve=curve,
        limitations=[
            "The equation is dR/dt = kin*C(t) - kout*R.",
            "C(t) is the accepted one-compartment IV bolus closed form.",
            "Multiplying kout by C(t) is not evaluated, because that integral is not elementary.",
            "The bolus concentration decays, so the steady state of this equation is 0.",
            "The evidence grade is predicted.",
            "Parameters were supplied by the caller and were not estimated.",
            "The curve is not a dosing recommendation and is not linked to a clinical endpoint.",
        ],
    )


def _response(r0: float, kin: float, kout: float, amount_over_volume: float, kel: float, time_h: float) -> float:
    decay = math.exp(-kout * time_h)
    gap = kout - kel
    if abs(gap) <= 1e-8 * max(kout, kel):
        driven = kin * amount_over_volume * time_h * decay
    else:
        driven = kin * amount_over_volume / gap * (math.exp(-kel * time_h) - decay)
    return r0 * decay + driven
