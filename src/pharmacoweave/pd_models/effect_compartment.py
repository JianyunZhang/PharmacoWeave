"""Effect-compartment concentration driven by the one-compartment IV bolus solution.

    dCe/dt = ke0 * (C - Ce)
    C(t) = (Dose / V) * exp(-(CL / V) * t)

Parameters are not estimated. The effect-site concentration stays predicted.
"""

import math
from typing import Literal

from pydantic import Field

from pharmacoweave.core.entities import Record
from pharmacoweave.errors import AnalysisError
from pharmacoweave.pk_models.analytical import iv_bolus_concentration


class EffectCompartmentPoint(Record):
    """Plasma and effect-site concentrations at one time. Both are predicted."""

    time_h: float
    plasma_concentration_mg_l: float
    effect_site_concentration_mg_l: float


class EffectCompartmentSimulation(Record):
    """A predicted effect-compartment curve for one IV bolus."""

    context_of_use: str = Field(min_length=1)
    evidence_grade: Literal["predicted"] = "predicted"
    dose_mg: float
    volume_l: float
    clearance_l_h: float
    ke0_per_h: float
    curve: list[EffectCompartmentPoint]
    limitations: list[str]
    constitutes_dosing_recommendation: Literal[False] = False
    schema_version: Literal["0.1.0"] = "0.1.0"


def simulate_effect_compartment(
    *,
    dose_mg: float,
    volume_l: float,
    clearance_l_h: float,
    ke0_per_h: float,
    times_h: list[float],
    context_of_use: str,
    estimate_parameters: bool = False,
) -> EffectCompartmentSimulation:
    """Evaluate the effect-compartment closed form. Estimation requests are rejected."""

    if estimate_parameters:
        raise AnalysisError("Effect-compartment parameters are inputs. This function does not estimate them.")
    if not math.isfinite(dose_mg) or dose_mg <= 0:
        raise AnalysisError("Dose must be positive.")
    if not math.isfinite(volume_l) or volume_l <= 0 or not math.isfinite(clearance_l_h) or clearance_l_h <= 0:
        raise AnalysisError("Volume and clearance must be positive.")
    if not math.isfinite(ke0_per_h) or ke0_per_h <= 0:
        raise AnalysisError("ke0 must be positive.")
    if len(times_h) < 1 or any(not math.isfinite(time) or time < 0 for time in times_h):
        raise AnalysisError("Effect-compartment times must be finite and zero or positive.")
    if any(later <= earlier for earlier, later in zip(times_h, times_h[1:], strict=False)):
        raise AnalysisError("Effect-compartment times must strictly increase.")
    kel = clearance_l_h / volume_l
    dose_over_volume = dose_mg / volume_l
    curve = [
        EffectCompartmentPoint(
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
            effect_site_concentration_mg_l=_effect_site(dose_over_volume, kel, ke0_per_h, time),
        )
        for time in times_h
    ]
    return EffectCompartmentSimulation(
        context_of_use=context_of_use,
        dose_mg=dose_mg,
        volume_l=volume_l,
        clearance_l_h=clearance_l_h,
        ke0_per_h=ke0_per_h,
        curve=curve,
        limitations=[
            "Plasma concentration is the accepted one-compartment IV bolus closed form.",
            "Effect-site concentration solves dCe/dt = ke0*(C-Ce) in closed form.",
            "ke0 equal to kel uses Ce(t) = (Dose/V)*kel*t*exp(-kel*t).",
            "The evidence grade is predicted. Effect-site concentration is not an observed concentration.",
            "Parameters were supplied by the caller and were not estimated.",
            "The curve is not a dosing recommendation.",
        ],
    )


def _effect_site(dose_over_volume: float, kel: float, ke0: float, time_h: float) -> float:
    if time_h == 0.0:
        return 0.0
    gap = ke0 - kel
    if abs(gap) <= 1e-8 * max(ke0, kel):
        return dose_over_volume * ke0 * time_h * math.exp(-kel * time_h)
    return dose_over_volume * (ke0 / gap) * (math.exp(-kel * time_h) - math.exp(-ke0 * time_h))
