"""Closed unit vocabulary and explicit molecular-weight conversions."""

import pytest

from pharmacoweave.core.entities import Compound
from pharmacoweave.core.units import (
    AUC_MASS,
    CLEARANCE,
    CONC_MASS,
    DOSE_PER_BODY_MASS,
    MASS_PER_AMOUNT,
    RATE,
    convert,
    convert_with_molecular_weight,
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


def test_dalton_and_kilodalton_match_mass_per_amount() -> None:
    assert parse_unit("Da").dimensions == MASS_PER_AMOUNT
    assert parse_unit("kDa").dimensions == MASS_PER_AMOUNT
    assert convert(1.0, "Da", "g/mol") == pytest.approx(1.0)
    assert convert(150.0, "kDa", "g/mol") == pytest.approx(150_000.0)


def test_explicit_molecular_weight_converts_mass_and_molar_quantities() -> None:
    # MW = 250 g/mol => 1 mg/L = 4 umol/L = 4000 nM
    assert convert(1.0, "mg/L", "mg/L") == pytest.approx(1.0)
    with pytest.raises(UnitError):
        convert(1.0, "mg/L", "uM")
    with pytest.raises(UnitError, match="molecular_weight"):
        convert_with_molecular_weight(1.0, "mg/L", "uM", molecular_weight=None)

    assert convert_with_molecular_weight(
        1.0, "mg/L", "uM", molecular_weight=250.0, molecular_weight_unit="g/mol"
    ) == pytest.approx(4.0)
    assert convert_with_molecular_weight(
        4000.0, "nM", "ug/mL", molecular_weight=250.0, molecular_weight_unit="Da"
    ) == pytest.approx(1.0)
    assert convert_with_molecular_weight(
        2.0, "mg/kg", "umol/kg", molecular_weight=250.0, molecular_weight_unit="g/mol"
    ) == pytest.approx(8.0)
    assert convert_with_molecular_weight(
        5.0, "mg*h/L", "umol*h/L", molecular_weight=250.0, molecular_weight_unit="g/mol"
    ) == pytest.approx(20.0)

    mab = Compound(compound_id="MAB-01", name="IgG1 mAb", molecular_weight=150.0, molecular_weight_unit="kDa")
    # 150 mg/L / 150,000 g/mol = 1e-6 mol/L = 1 uM
    assert mab.convert_quantity(150.0, "mg/L", "uM") == pytest.approx(1.0)
    assert mab.convert_quantity(1.0, "uM", "ug/mL") == pytest.approx(150.0)

    no_mw = Compound(compound_id="CMP-NO-MW", name="Unweighed")
    with pytest.raises(UnitError, match="molecular_weight"):
        no_mw.convert_quantity(1.0, "uM", "mg/L")
