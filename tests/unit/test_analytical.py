"""Analytical PK reference values across 1-compartment and 2-compartment models."""

import math

import pytest

from pharmacoweave.errors import UnitError
from pharmacoweave.pk_models import (
    BUILTIN_SPECS,
    ONE_COMPARTMENT_EXTRAVASCULAR,
    ONE_COMPARTMENT_IV_BOLUS,
    TWO_COMPARTMENT_IV_BOLUS,
    auc_inf_iv_bolus,
    half_life,
    iv_bolus_concentration,
    iv_infusion_concentration,
    one_compartment_oral_concentration,
    two_compartment_iv_bolus_concentration,
)


def test_bolus_concentration_and_auc() -> None:
    concentration = iv_bolus_concentration(
        dose=0.25,
        dose_unit="mg",
        volume=0.25,
        volume_unit="L",
        clearance=0.05,
        clearance_unit="L/h",
        time=0.0,
        time_unit="h",
        concentration_unit="ug/mL",
    )
    assert concentration == pytest.approx(1.0)
    assert half_life(volume=0.25, volume_unit="L", clearance=0.05, clearance_unit="L/h") == pytest.approx(
        math.log(2.0) / 0.2
    )
    assert auc_inf_iv_bolus(
        dose=0.25,
        dose_unit="mg",
        clearance=0.05,
        clearance_unit="L/h",
        auc_unit="ug*h/mL",
    ) == pytest.approx(5.0)


def test_clearance_can_be_supplied_in_milliliters_per_minute() -> None:
    concentration = iv_bolus_concentration(
        dose=0.25,
        dose_unit="mg",
        volume=250.0,
        volume_unit="mL",
        clearance=0.05 / 0.06,
        clearance_unit="mL/min",
        time=60.0,
        time_unit="min",
        concentration_unit="ng/mL",
    )
    assert concentration == pytest.approx(1000.0 * math.exp(-0.2))


def test_reference_rejects_incompatible_units_and_negative_time() -> None:
    with pytest.raises(UnitError):
        iv_bolus_concentration(
            dose=1.0,
            dose_unit="h",
            volume=1.0,
            volume_unit="L",
            clearance=1.0,
            clearance_unit="L/h",
            time=0.0,
            time_unit="h",
            concentration_unit="mg/L",
        )
    with pytest.raises(ValueError):
        iv_bolus_concentration(
            dose=1.0,
            dose_unit="mg",
            volume=1.0,
            volume_unit="L",
            clearance=1.0,
            clearance_unit="L/h",
            time=-1.0,
            time_unit="h",
            concentration_unit="mg/L",
        )


def test_infusion_and_bateman_and_two_compartment_analytical_solutions() -> None:
    c_end = iv_infusion_concentration(
        dose=2.0,
        dose_unit="mg",
        infusion_duration=1.0,
        infusion_duration_unit="h",
        volume=1.0,
        volume_unit="L",
        clearance=0.25,
        clearance_unit="L/h",
        time=1.0,
        time_unit="h",
        concentration_unit="mg/L",
    )
    assert c_end == pytest.approx((2.0 / 0.25) * (1.0 - math.exp(-0.25)))
    c_post = iv_infusion_concentration(
        dose=2.0,
        dose_unit="mg",
        infusion_duration=1.0,
        infusion_duration_unit="h",
        volume=1.0,
        volume_unit="L",
        clearance=0.25,
        clearance_unit="L/h",
        time=3.0,
        time_unit="h",
        concentration_unit="mg/L",
    )
    assert c_post == pytest.approx(c_end * math.exp(-0.5))

    # Confluent ka == ke limit: C(t) = (Dose / (V/F)) * ka * t * exp(-ke * t)
    c_confluent = one_compartment_oral_concentration(
        dose=2.0,
        dose_unit="mg",
        ka=0.5,
        ka_unit="1/h",
        apparent_volume=1.0,
        volume_unit="L",
        apparent_clearance=0.5,
        clearance_unit="L/h",
        time=2.0,
        time_unit="h",
        concentration_unit="uM",
        molecular_weight=250.0,
        molecular_weight_unit="g/mol",
    )
    assert c_confluent == pytest.approx((2.0 * 0.5 * 2.0 * math.exp(-1.0)) * 4.0)

    c0_2c = two_compartment_iv_bolus_concentration(
        dose=10.0,
        dose_unit="mg",
        clearance=0.5,
        clearance_unit="L/h",
        central_volume=2.0,
        central_volume_unit="L",
        intercompartmental_clearance=1.5,
        intercompartmental_clearance_unit="L/h",
        peripheral_volume=4.0,
        peripheral_volume_unit="L",
        time=0.0,
        time_unit="h",
        concentration_unit="mg/L",
    )
    assert c0_2c == pytest.approx(5.0)


def test_builtin_specs_cover_all_implemented_compartmental_models() -> None:
    assert ONE_COMPARTMENT_IV_BOLUS.estimation_implemented is True
    assert ONE_COMPARTMENT_IV_BOLUS.parameters == ["CL", "V"]
    assert "Dose/V" in ONE_COMPARTMENT_IV_BOLUS.equations_ref
    assert "1/y^2" in " ".join(ONE_COMPARTMENT_IV_BOLUS.assumptions)

    assert ONE_COMPARTMENT_EXTRAVASCULAR.estimation_implemented is True
    assert ONE_COMPARTMENT_EXTRAVASCULAR.parameters == ["ka", "CL/F", "V/F"]

    assert TWO_COMPARTMENT_IV_BOLUS.estimation_implemented is True
    assert TWO_COMPARTMENT_IV_BOLUS.parameters == ["CL", "V1", "Q", "V2"]
    assert set(BUILTIN_SPECS) == {
        "one_compartment_iv_bolus",
        "one_compartment_extravascular",
        "two_compartment_iv_bolus",
    }
