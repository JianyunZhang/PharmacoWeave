"""Analytical IV bolus reference values."""

import math

import pytest

from pharmacoweave.errors import UnitError
from pharmacoweave.pk_models import ONE_COMPARTMENT_IV_BOLUS, auc_inf_iv_bolus, half_life, iv_bolus_concentration


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
    # 0.05 L/h = 0.05 / 0.06 mL/min
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


def test_builtin_spec_is_not_an_estimator() -> None:
    assert ONE_COMPARTMENT_IV_BOLUS.estimation_implemented is False
    assert ONE_COMPARTMENT_IV_BOLUS.parameters == ["CL", "V"]
    assert "Dose/V" in ONE_COMPARTMENT_IV_BOLUS.equations_ref
