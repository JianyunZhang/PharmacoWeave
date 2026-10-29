"""Analytical one-compartment IV bolus reference.

These functions generate and check concentrations from known parameters.
Estimation of CL and V lives in ``pharmacoweave.pk_models.fit``.
"""

import math

from pharmacoweave.core.units import convert


def iv_bolus_concentration(
    *,
    dose: float,
    dose_unit: str,
    volume: float,
    volume_unit: str,
    clearance: float,
    clearance_unit: str,
    time: float,
    time_unit: str,
    concentration_unit: str,
) -> float:
    """Return C(t) = (Dose / V) * exp(-(CL / V) * t).

    Mass and molar concentrations are not interconverted. Pass a molecular
    weight through a separate, explicit step before requesting molar units.
    """

    if dose < 0 or volume <= 0 or clearance <= 0 or time < 0:
        raise ValueError("dose must be >= 0, volume and clearance must be > 0, and time must be >= 0.")
    dose_mg = convert(dose, dose_unit, "mg")
    volume_l = convert(volume, volume_unit, "L")
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    time_h = convert(time, time_unit, "h")
    elimination_rate = clearance_l_h / volume_l
    concentration_mg_l = (dose_mg / volume_l) * math.exp(-elimination_rate * time_h)
    return convert(concentration_mg_l, "mg/L", concentration_unit)


def half_life(
    *,
    volume: float,
    volume_unit: str,
    clearance: float,
    clearance_unit: str,
    time_unit: str = "h",
) -> float:
    """Return ln(2) * V / CL for a one-compartment linear model."""

    if volume <= 0 or clearance <= 0:
        raise ValueError("volume and clearance must be positive.")
    volume_l = convert(volume, volume_unit, "L")
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    half_life_h = math.log(2.0) * volume_l / clearance_l_h
    return convert(half_life_h, "h", time_unit)


def auc_inf_iv_bolus(
    *,
    dose: float,
    dose_unit: str,
    clearance: float,
    clearance_unit: str,
    auc_unit: str,
) -> float:
    """Return Dose / CL, the IV bolus AUC from time zero to infinity."""

    if dose < 0 or clearance <= 0:
        raise ValueError("dose must be >= 0 and clearance must be positive.")
    dose_mg = convert(dose, dose_unit, "mg")
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    auc_mg_h_l = dose_mg / clearance_l_h
    return convert(auc_mg_h_l, "mg*h/L", auc_unit)
