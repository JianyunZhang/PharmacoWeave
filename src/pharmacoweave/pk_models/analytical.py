"""Analytical compartmental PK reference solutions.

These functions generate and check concentrations from known parameters.
Mass and molar units stay separate unless ``molecular_weight`` is explicitly
supplied.
"""

import math

from pharmacoweave.core.units import convert, convert_with_molecular_weight


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
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return C(t) = (Dose / V) * exp(-(CL / V) * t)."""

    if dose < 0 or volume <= 0 or clearance <= 0 or time < 0:
        raise ValueError("dose must be >= 0, volume and clearance must be > 0, and time must be >= 0.")
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    volume_l = convert(volume, volume_unit, "L")
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    time_h = convert(time, time_unit, "h")
    elimination_rate = clearance_l_h / volume_l
    concentration_mg_l = (dose_mg / volume_l) * math.exp(-elimination_rate * time_h)
    return convert_with_molecular_weight(
        concentration_mg_l,
        "mg/L",
        concentration_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def iv_infusion_concentration(
    *,
    dose: float,
    dose_unit: str,
    infusion_duration: float,
    infusion_duration_unit: str,
    volume: float,
    volume_unit: str,
    clearance: float,
    clearance_unit: str,
    time: float,
    time_unit: str,
    concentration_unit: str,
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return one-compartment constant-rate IV infusion concentration C(t)."""

    if dose < 0 or infusion_duration <= 0 or volume <= 0 or clearance <= 0 or time < 0:
        raise ValueError(
            "dose must be >= 0, infusion_duration, volume, and clearance must be > 0, and time must be >= 0."
        )
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    duration_h = convert(infusion_duration, infusion_duration_unit, "h")
    volume_l = convert(volume, volume_unit, "L")
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    time_h = convert(time, time_unit, "h")
    if dose_mg == 0 or time_h == 0:
        return convert_with_molecular_weight(
            0.0,
            "mg/L",
            concentration_unit,
            molecular_weight=molecular_weight,
            molecular_weight_unit=molecular_weight_unit,
        )
    rate_mg_h = dose_mg / duration_h
    kel = clearance_l_h / volume_l
    active_h = min(time_h, duration_h)
    at_end_or_now = -(rate_mg_h / clearance_l_h) * math.expm1(-kel * active_h)
    concentration_mg_l = (
        at_end_or_now if time_h <= duration_h else at_end_or_now * math.exp(-kel * (time_h - duration_h))
    )
    return convert_with_molecular_weight(
        concentration_mg_l,
        "mg/L",
        concentration_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def one_compartment_oral_concentration(
    *,
    dose: float,
    dose_unit: str,
    ka: float,
    ka_unit: str,
    apparent_volume: float,
    volume_unit: str,
    apparent_clearance: float,
    clearance_unit: str,
    time: float,
    time_unit: str,
    concentration_unit: str,
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return Bateman extravascular C(t) parameterized by ka, CL/F, and V/F."""

    if dose < 0 or ka <= 0 or apparent_volume <= 0 or apparent_clearance <= 0 or time < 0:
        raise ValueError(
            "dose must be >= 0, ka, apparent_volume, and apparent_clearance must be > 0, and time must be >= 0."
        )
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    ka_h = convert(ka, ka_unit, "1/h")
    vf_l = convert(apparent_volume, volume_unit, "L")
    clf_l_h = convert(apparent_clearance, clearance_unit, "L/h")
    time_h = convert(time, time_unit, "h")
    concentration_mg_l = bateman_mg_l(dose_mg=dose_mg, ka_h=ka_h, clf_l_h=clf_l_h, vf_l=vf_l, time_h=time_h)
    return convert_with_molecular_weight(
        concentration_mg_l,
        "mg/L",
        concentration_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def two_compartment_iv_bolus_concentration(
    *,
    dose: float,
    dose_unit: str,
    clearance: float,
    clearance_unit: str,
    central_volume: float,
    central_volume_unit: str,
    intercompartmental_clearance: float,
    intercompartmental_clearance_unit: str,
    peripheral_volume: float,
    peripheral_volume_unit: str,
    time: float,
    time_unit: str,
    concentration_unit: str,
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return central-compartment C(t) = A*exp(-alpha*t) + B*exp(-beta*t) for CL, V1, Q, V2."""

    if (
        dose < 0
        or clearance <= 0
        or central_volume <= 0
        or intercompartmental_clearance <= 0
        or peripheral_volume <= 0
        or time < 0
    ):
        raise ValueError(
            "dose must be >= 0, clearance, central_volume, intercompartmental_clearance, "
            "and peripheral_volume must be > 0, and time must be >= 0."
        )
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    cl_l_h = convert(clearance, clearance_unit, "L/h")
    v1_l = convert(central_volume, central_volume_unit, "L")
    q_l_h = convert(intercompartmental_clearance, intercompartmental_clearance_unit, "L/h")
    v2_l = convert(peripheral_volume, peripheral_volume_unit, "L")
    time_h = convert(time, time_unit, "h")
    concentration_mg_l = two_compartment_bolus_mg_l(
        dose_mg=dose_mg,
        cl_l_h=cl_l_h,
        v1_l=v1_l,
        q_l_h=q_l_h,
        v2_l=v2_l,
        time_h=time_h,
    )
    return convert_with_molecular_weight(
        concentration_mg_l,
        "mg/L",
        concentration_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def bateman_mg_l(*, dose_mg: float, ka_h: float, clf_l_h: float, vf_l: float, time_h: float) -> float:
    """Evaluate Bateman C(t) in canonical mg/L without cancellation when ka ≈ ke."""

    if dose_mg == 0.0 or time_h == 0.0:
        return 0.0
    kel = clf_l_h / vf_l
    diff = ka_h - kel
    delta = diff * time_h
    exp_kel = math.exp(-kel * time_h)
    if abs(delta) <= 1e-7:
        series = time_h * (1.0 - 0.5 * delta + (delta * delta) / 6.0 - (delta * delta * delta) / 24.0)
        return (dose_mg / vf_l) * ka_h * exp_kel * series
    return (dose_mg / vf_l) * (ka_h / diff) * exp_kel * (-math.expm1(-delta))


def two_compartment_macro_constants(
    *,
    dose_mg: float,
    cl_l_h: float,
    v1_l: float,
    q_l_h: float,
    v2_l: float,
) -> tuple[float, float, float, float]:
    """Return (A_mg_l, alpha_h, B_mg_l, beta_h) with alpha > beta > 0."""

    k10 = cl_l_h / v1_l
    k12 = q_l_h / v1_l
    k21 = q_l_h / v2_l
    rate_sum = k10 + k12 + k21
    rate_prod = k10 * k21
    disc_inner = (k10 + k12 - k21) * (k10 + k12 - k21) + 4.0 * k12 * k21
    root = math.sqrt(max(0.0, disc_inner))
    alpha = 0.5 * (rate_sum + root)
    beta = (2.0 * rate_prod) / (rate_sum + root)
    c0 = dose_mg / v1_l
    span = alpha - beta
    a_coef = c0 * (alpha - k21) / span
    b_coef = c0 * (k21 - beta) / span
    return a_coef, alpha, b_coef, beta


def two_compartment_bolus_mg_l(
    *,
    dose_mg: float,
    cl_l_h: float,
    v1_l: float,
    q_l_h: float,
    v2_l: float,
    time_h: float,
) -> float:
    """Evaluate two-compartment IV bolus C(t) in canonical mg/L."""

    if dose_mg == 0.0:
        return 0.0
    a_coef, alpha, b_coef, beta = two_compartment_macro_constants(
        dose_mg=dose_mg,
        cl_l_h=cl_l_h,
        v1_l=v1_l,
        q_l_h=q_l_h,
        v2_l=v2_l,
    )
    return a_coef * math.exp(-alpha * time_h) + b_coef * math.exp(-beta * time_h)


def three_compartment_iv_bolus_concentration(
    *,
    dose: float,
    dose_unit: str,
    clearance: float,
    clearance_unit: str,
    central_volume: float,
    central_volume_unit: str,
    q2: float,
    q2_unit: str,
    peripheral_volume_2: float,
    peripheral_volume_2_unit: str,
    q3: float,
    q3_unit: str,
    peripheral_volume_3: float,
    peripheral_volume_3_unit: str,
    time: float,
    time_unit: str,
    concentration_unit: str,
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return central C(t) for CL, V1, Q2, V2, Q3, and V3."""

    if (
        dose < 0
        or clearance <= 0
        or central_volume <= 0
        or q2 <= 0
        or peripheral_volume_2 <= 0
        or q3 <= 0
        or peripheral_volume_3 <= 0
        or time < 0
    ):
        raise ValueError(
            "dose must be >= 0, clearance, volumes, and intercompartmental clearances must be > 0, "
            "and time must be >= 0."
        )
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    concentration_mg_l = three_compartment_bolus_mg_l(
        dose_mg=dose_mg,
        cl_l_h=convert(clearance, clearance_unit, "L/h"),
        v1_l=convert(central_volume, central_volume_unit, "L"),
        q2_l_h=convert(q2, q2_unit, "L/h"),
        v2_l=convert(peripheral_volume_2, peripheral_volume_2_unit, "L"),
        q3_l_h=convert(q3, q3_unit, "L/h"),
        v3_l=convert(peripheral_volume_3, peripheral_volume_3_unit, "L"),
        time_h=convert(time, time_unit, "h"),
    )
    return convert_with_molecular_weight(
        concentration_mg_l,
        "mg/L",
        concentration_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )


def three_compartment_macro_constants(
    *,
    dose_mg: float,
    cl_l_h: float,
    v1_l: float,
    q2_l_h: float,
    v2_l: float,
    q3_l_h: float,
    v3_l: float,
) -> tuple[float, float, float, float, float, float]:
    """Return (A, alpha, B, beta, G, gamma) in mg/L and 1/h, with alpha > beta > gamma > 0."""

    k10 = cl_l_h / v1_l
    k12 = q2_l_h / v1_l
    k21 = q2_l_h / v2_l
    k13 = q3_l_h / v1_l
    k31 = q3_l_h / v3_l
    alpha, beta, gamma = _three_decay_rates(k10, k12, k21, k13, k31)
    c0 = dose_mg / v1_l
    a_coef = c0 * _hybrid_fraction(alpha, beta, gamma, k21, k31)
    b_coef = c0 * _hybrid_fraction(beta, alpha, gamma, k21, k31)
    g_coef = c0 * _hybrid_fraction(gamma, alpha, beta, k21, k31)
    return a_coef, alpha, b_coef, beta, g_coef, gamma


def three_compartment_bolus_mg_l(
    *,
    dose_mg: float,
    cl_l_h: float,
    v1_l: float,
    q2_l_h: float,
    v2_l: float,
    q3_l_h: float,
    v3_l: float,
    time_h: float,
) -> float:
    """Evaluate three-compartment IV bolus C(t) in canonical mg/L."""

    if dose_mg == 0.0:
        return 0.0
    a_coef, alpha, b_coef, beta, g_coef, gamma = three_compartment_macro_constants(
        dose_mg=dose_mg,
        cl_l_h=cl_l_h,
        v1_l=v1_l,
        q2_l_h=q2_l_h,
        v2_l=v2_l,
        q3_l_h=q3_l_h,
        v3_l=v3_l,
    )
    return (
        a_coef * math.exp(-alpha * time_h)
        + b_coef * math.exp(-beta * time_h)
        + g_coef * math.exp(-gamma * time_h)
    )


def _three_decay_rates(k10: float, k12: float, k21: float, k13: float, k31: float) -> tuple[float, float, float]:
    elimination = k10 + k12 + k13
    sum_rates = elimination + k21 + k31
    sum_products = k10 * (k21 + k31) + k21 * k31 + k12 * k31 + k13 * k21
    product = k10 * k21 * k31
    roots = _cubic_real_roots(-sum_rates, sum_products, -product)
    positive = sorted((root for root in roots if root > 1e-12), reverse=True)
    if len(positive) != 3 or positive[0] - positive[1] < 1e-8 or positive[1] - positive[2] < 1e-8:
        raise ValueError("Three-compartment hybrid rates are not three separated positive values.")
    return positive[0], positive[1], positive[2]


def _hybrid_fraction(rate: float, other_a: float, other_b: float, k21: float, k31: float) -> float:
    return ((rate - k21) * (rate - k31)) / ((rate - other_a) * (rate - other_b))


def _cubic_real_roots(linear2: float, linear1: float, constant: float) -> tuple[float, float, float]:
    """Return three real roots of x^3 + linear2 x^2 + linear1 x + constant = 0."""

    shift = linear2 / 3.0
    depressed_a = linear1 - linear2 * linear2 / 3.0
    depressed_b = (2.0 * linear2**3 - 9.0 * linear2 * linear1 + 27.0 * constant) / 27.0
    if depressed_a >= 0.0:
        raise ValueError("Three-compartment characteristic polynomial does not have three real roots.")
    rho = math.sqrt(-depressed_a / 3.0)
    cosine = (-depressed_b / 2.0) / (rho**3)
    cosine = max(-1.0, min(1.0, cosine))
    theta = math.acos(cosine)
    first = 2.0 * rho * math.cos(theta / 3.0)
    second = 2.0 * rho * math.cos((theta + 2.0 * math.pi) / 3.0)
    third = 2.0 * rho * math.cos((theta + 4.0 * math.pi) / 3.0)
    return (first - shift, second - shift, third - shift)


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
    molecular_weight: float | None = None,
    molecular_weight_unit: str | None = "g/mol",
) -> float:
    """Return Dose / CL, the IV AUC from time zero to infinity."""

    if dose < 0 or clearance <= 0:
        raise ValueError("dose must be >= 0 and clearance must be positive.")
    dose_mg = convert_with_molecular_weight(
        dose,
        dose_unit,
        "mg",
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
    clearance_l_h = convert(clearance, clearance_unit, "L/h")
    auc_mg_h_l = dose_mg / clearance_l_h
    return convert_with_molecular_weight(
        auc_mg_h_l,
        "mg*h/L",
        auc_unit,
        molecular_weight=molecular_weight,
        molecular_weight_unit=molecular_weight_unit,
    )
