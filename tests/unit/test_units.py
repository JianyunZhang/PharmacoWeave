"""Closed unit vocabulary."""

import pytest

from pharmacoweave.core.units import (
    AUC_MASS,
    CLEARANCE,
    CONC_MASS,
    DOSE_PER_BODY_MASS,
    RATE,
    convert,
    parse_unit,
)
from pharmacoweave.errors import UnitError


def test_microgram_per_milliliter_matches_milligram_per_liter() -> None:
    assert convert(1.0, "ug/mL", "mg/L") == pytest.approx(1.0)
    assert convert(1.0, "µg/mL", "mg/L") == pytest.approx(1.0)
    assert convert(1.0, "mcg/mL", "mg/L") == pytest.approx(1.0)


def test_clearance_minute_to_hour() -> None:
    assert convert(1.0, "mL/min", "L/h") == pytest.approx(0.06)


def test_body_weight_dose_is_not_a_drug_mass() -> None:
    assert parse_unit("mg/kg").dimensions == DOSE_PER_BODY_MASS
    with pytest.raises(UnitError):
        convert(1.0, "mg/kg", "mg")


def test_kilogram_is_body_mass_and_does_not_convert_to_drug_grams() -> None:
    with pytest.raises(UnitError):
        convert(1.0, "kg", "g")


def test_auc_and_rate_dimensions() -> None:
    assert parse_unit("mg*h/L").dimensions == AUC_MASS
    assert parse_unit("1/h").dimensions == RATE
    assert parse_unit("ug*h/mL").dimensions == AUC_MASS
    assert convert(5.0, "mg*h/L", "ug*h/mL") == pytest.approx(5.0)


def test_infusion_rate_keeps_body_mass_distinct() -> None:
    parsed = parse_unit("mg/kg/h")
    assert parsed.dimensions == (-1, 1, -1, 0, 0)
    assert parsed.dimensions != CLEARANCE
    assert parsed.dimensions != CONC_MASS


def test_unknown_unit_is_rejected() -> None:
    with pytest.raises(UnitError):
        parse_unit("foot")


def test_percent_is_dimensionless_and_is_not_a_concentration() -> None:
    assert convert(50.0, "percent", "1") == pytest.approx(0.5)
    assert convert(1.0, "fold", "1") == pytest.approx(1.0)
    with pytest.raises(UnitError):
        convert(50.0, "percent", "mg/L")
