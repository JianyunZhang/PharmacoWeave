"""Closed unit vocabulary for PK quantities.

Drug mass and body mass are different dimensions, so ``mg/kg`` does not
collapse to a number. The atom ``kg`` means body mass. Drug amounts use
``pg``, ``ng``, ``ug`` / ``µg`` / ``mcg``, ``mg``, or ``g``.

Canonical scale factors are: hour, milligram of drug, kilogram of body
weight, litre, and mole. Conversions inside one dimension are linear.
Molar and mass concentrations are not converted into each other, because
that requires a molecular weight and has to be an explicit step.

``percent`` and ``pct`` are a dimensionless effect scale (50 percent = 0.5).
``fold`` is a dimensionless factor of 1, the same dimension as ``1``.
Neither is a concentration.
"""

import math
import re
from dataclasses import dataclass
from typing import TypeAlias

from pharmacoweave.errors import UnitError

Dimensions: TypeAlias = tuple[int, int, int, int, int]

TIME: Dimensions = (1, 0, 0, 0, 0)
DRUG_MASS: Dimensions = (0, 1, 0, 0, 0)
BODY_MASS: Dimensions = (0, 0, 1, 0, 0)
VOLUME: Dimensions = (0, 0, 0, 1, 0)
AMOUNT: Dimensions = (0, 0, 0, 0, 1)
CONC_MASS: Dimensions = (0, 1, 0, -1, 0)
CONC_AMOUNT: Dimensions = (0, 0, 0, -1, 1)
CLEARANCE: Dimensions = (-1, 0, 0, 1, 0)
DOSE_PER_BODY_MASS: Dimensions = (0, 1, -1, 0, 0)
AMOUNT_PER_BODY_MASS: Dimensions = (0, 0, -1, 0, 1)
AUC_MASS: Dimensions = (1, 1, 0, -1, 0)
AUC_AMOUNT: Dimensions = (1, 0, 0, -1, 1)
RATE: Dimensions = (-1, 0, 0, 0, 0)
MASS_PER_AMOUNT: Dimensions = (0, 1, 0, 0, -1)
DIMENSIONLESS: Dimensions = (0, 0, 0, 0, 0)

_NAMES = ("time", "drug_mass", "body_mass", "volume", "amount")

# factor converts one atom into the canonical unit of its dimension.
_ATOMS: dict[str, tuple[Dimensions, float]] = {
    "s": (TIME, 1.0 / 3600.0),
    "sec": (TIME, 1.0 / 3600.0),
    "second": (TIME, 1.0 / 3600.0),
    "seconds": (TIME, 1.0 / 3600.0),
    "min": (TIME, 1.0 / 60.0),
    "mins": (TIME, 1.0 / 60.0),
    "minute": (TIME, 1.0 / 60.0),
    "minutes": (TIME, 1.0 / 60.0),
    "h": (TIME, 1.0),
    "hr": (TIME, 1.0),
    "hrs": (TIME, 1.0),
    "hour": (TIME, 1.0),
    "hours": (TIME, 1.0),
    "d": (TIME, 24.0),
    "day": (TIME, 24.0),
    "days": (TIME, 24.0),
    "pg": (DRUG_MASS, 1e-9),
    "ng": (DRUG_MASS, 1e-6),
    "ug": (DRUG_MASS, 1e-3),
    "mcg": (DRUG_MASS, 1e-3),
    "mg": (DRUG_MASS, 1.0),
    "g": (DRUG_MASS, 1e3),
    "kg": (BODY_MASS, 1.0),
    "uL": (VOLUME, 1e-6),
    "mL": (VOLUME, 1e-3),
    "L": (VOLUME, 1.0),
    "liter": (VOLUME, 1.0),
    "litre": (VOLUME, 1.0),
    "liters": (VOLUME, 1.0),
    "litres": (VOLUME, 1.0),
    "pmol": (AMOUNT, 1e-12),
    "nmol": (AMOUNT, 1e-9),
    "umol": (AMOUNT, 1e-6),
    "mmol": (AMOUNT, 1e-3),
    "mol": (AMOUNT, 1.0),
    "1": (DIMENSIONLESS, 1.0),
    "fold": (DIMENSIONLESS, 1.0),
    "percent": (DIMENSIONLESS, 0.01),
    "pct": (DIMENSIONLESS, 0.01),
}

_TOKEN = re.compile(r"[A-Za-z]+|1")


@dataclass(frozen=True)
class Unit:
    """A parsed unit expression."""

    expression: str
    dimensions: Dimensions
    factor: float


def parse_unit(expression: str) -> Unit:
    """Parse a unit expression such as ``ug/mL`` or ``mg*h/L``."""

    if not isinstance(expression, str) or not expression.strip():
        raise UnitError("A unit expression is required.")
    text = expression.strip().replace("µ", "u").replace("μ", "u").replace("·", "*").replace(" ", "")
    tokens = _tokenize(text, expression)
    dimensions = [0, 0, 0, 0, 0]
    factor = 1.0
    expect_atom = True
    operator = "*"
    for token in tokens:
        if expect_atom:
            atom_dimensions, atom_factor = _lookup(token, expression)
            if operator == "*":
                factor *= atom_factor
                dimensions = [left + right for left, right in zip(dimensions, atom_dimensions, strict=True)]
            else:
                factor /= atom_factor
                dimensions = [left - right for left, right in zip(dimensions, atom_dimensions, strict=True)]
            expect_atom = False
            continue
        if token not in {"*", "/"}:
            raise UnitError(f"Expected * or / in unit {expression!r}, found {token!r}.")
        operator = token
        expect_atom = True
    if expect_atom:
        raise UnitError(f"Unit {expression!r} ends with an operator.")
    parsed: Dimensions = (
        dimensions[0],
        dimensions[1],
        dimensions[2],
        dimensions[3],
        dimensions[4],
    )
    return Unit(expression=expression.strip(), dimensions=parsed, factor=factor)


def convert(value: float, source: str, target: str) -> float:
    """Convert ``value`` from ``source`` to ``target`` inside one dimension."""

    if not math.isfinite(value):
        raise UnitError("Cannot convert a non-finite value.")
    parsed_source = parse_unit(source)
    parsed_target = parse_unit(target)
    if parsed_source.dimensions != parsed_target.dimensions:
        raise UnitError(
            f"Cannot convert {source} ({describe(parsed_source.dimensions)}) "
            f"to {target} ({describe(parsed_target.dimensions)})."
        )
    return value * parsed_source.factor / parsed_target.factor


def describe(dimensions: Dimensions) -> str:
    """Return a readable dimension label for error messages."""

    numerator: list[str] = []
    denominator: list[str] = []
    for name, power in zip(_NAMES, dimensions, strict=True):
        if power > 0:
            numerator.append(name if power == 1 else f"{name}^{power}")
        elif power < 0:
            magnitude = -power
            denominator.append(name if magnitude == 1 else f"{name}^{magnitude}")
    if not numerator and not denominator:
        return "dimensionless"
    left = "*".join(numerator) if numerator else "1"
    if not denominator:
        return left
    return f"{left}/{'*'.join(denominator)}"


def same_dimensions(unit: str, expected: Dimensions) -> bool:
    """Return whether ``unit`` has the expected dimension."""

    return parse_unit(unit).dimensions == expected


def _tokenize(text: str, original: str) -> list[str]:
    tokens: list[str] = []
    index = 0
    while index < len(text):
        character = text[index]
        if character in {"*", "/"}:
            tokens.append(character)
            index += 1
            continue
        match = _TOKEN.match(text, index)
        if match is None:
            raise UnitError(f"Cannot parse unit {original!r} at {text[index:]!r}.")
        tokens.append(match.group(0))
        index = match.end()
    if not tokens:
        raise UnitError(f"Cannot parse unit {original!r}.")
    return tokens


def _lookup(token: str, original: str) -> tuple[Dimensions, float]:
    try:
        return _ATOMS[token]
    except KeyError as exc:
        raise UnitError(
            f"Unknown unit atom {token!r} in {original!r}. Use the PharmacoWeave unit vocabulary; kg means body mass."
        ) from exc
