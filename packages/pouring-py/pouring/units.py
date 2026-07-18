"""
Units and dimensional checking (spec 5.2-5.3).

Decimal throughout, never binary float. Spec 5.1: floating point breaks
equality, hashing, canonical serialisation, and significant figures, and a
measurement written as "0.100" means something different from "0.1".
Simulation converts to float at its own boundary; quantities do not.

Two conversions are traps, and both are handled explicitly rather than by the
usual multiply-by-a-factor:

  * **Celsius is affine.** 20 degC is 293.15 K, not 20 K scaled by anything.
  * **Gauge pressure is offset from absolute.** 1 barg is 2.013 bar absolute.
    Treating them as the same number is how vessels get over-pressurised.

Fractions are deliberately *not* interconvertible. Mass fraction, mole
fraction, and volume fraction are dimensionless but not the same quantity, and
converting between them needs composition data this module does not have. They
are therefore separate dimensions, so the attempt fails loudly.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from typing import Optional

# Standard atmosphere, for gauge-to-absolute pressure. An assumption, and one
# worth naming: gauge readings are relative to ambient, which is not always
# exactly this.
STANDARD_ATMOSPHERE = Decimal("101325")


@dataclass(frozen=True)
class UnitDefinition:
    symbol: str
    dimension: str
    # base = value * scale + offset
    scale: Decimal
    offset: Decimal = Decimal(0)


def _u(symbol: str, dimension: str, scale: str, offset: str = "0") -> UnitDefinition:
    return UnitDefinition(symbol, dimension, Decimal(scale), Decimal(offset))


UNITS: dict[str, UnitDefinition] = {
    # amount — base mol
    "mol": _u("mol", "amount", "1"),
    "mmol": _u("mmol", "amount", "0.001"),
    "umol": _u("umol", "amount", "0.000001"),
    # mass — base g
    "g": _u("g", "mass", "1"),
    "mg": _u("mg", "mass", "0.001"),
    "kg": _u("kg", "mass", "1000"),
    # volume — base L
    "L": _u("L", "volume", "1"),
    "mL": _u("mL", "volume", "0.001"),
    "uL": _u("uL", "volume", "0.000001"),
    # temperature — base K, affine
    "K": _u("K", "temperature", "1"),
    "degC": _u("degC", "temperature", "1", "273.15"),
    # time — base s
    "s": _u("s", "time", "1"),
    "min": _u("min", "time", "60"),
    "h": _u("h", "time", "3600"),
    # concentration — base M
    "M": _u("M", "concentration", "1"),
    "mM": _u("mM", "concentration", "0.001"),
    "uM": _u("uM", "concentration", "0.000001"),
    "nM": _u("nM", "concentration", "0.000000001"),
    # pressure — base Pa absolute
    "Pa": _u("Pa", "pressure", "1"),
    "kPa": _u("kPa", "pressure", "1000"),
    "bar": _u("bar", "pressure", "100000"),
    "atm": _u("atm", "pressure", "101325"),
    "mmHg": _u("mmHg", "pressure", "133.322"),
    # pressure, gauge — offset from absolute by one atmosphere
    "barg": _u("barg", "pressure", "100000", str(STANDARD_ATMOSPHERE)),
    "psig": _u("psig", "pressure", "6894.757", str(STANDARD_ATMOSPHERE)),
    # flow rate — base L/s
    "L/s": _u("L/s", "flowRate", "1"),
    "mL/min": _u("mL/min", "flowRate", "0.0000166666666666667"),
    "L/h": _u("L/h", "flowRate", "0.000277777777777778"),
    # stirring — base rpm
    "rpm": _u("rpm", "stirringRate", "1"),
    # length — base m
    "m": _u("m", "length", "1"),
    "cm": _u("cm", "length", "0.01"),
    "mm": _u("mm", "length", "0.001"),
    # density — base g/mL
    "g/mL": _u("g/mL", "density", "1"),
    "kg/L": _u("kg/L", "density", "1"),
    # fractions — dimensionless but NOT interchangeable
    "massFraction": _u("massFraction", "massFraction", "1"),
    "moleFraction": _u("moleFraction", "moleFraction", "1"),
    "volumeFraction": _u("volumeFraction", "volumeFraction", "1"),
    # misc
    "equiv": _u("equiv", "equivalents", "1"),
    "nm": _u("nm", "wavelength", "1"),
    "V": _u("V", "potential", "1"),
    "A": _u("A", "current", "1"),
    "pH": _u("pH", "acidity", "1"),
}


class UnknownUnit(Exception):
    """Raised for a unit not in the registry — diagnostic `unit.unknown`."""


class DimensionMismatch(Exception):
    """Raised across dimensions — diagnostic `unit.dimension`."""


def is_known_unit(symbol: str) -> bool:
    return symbol in UNITS


def definition(symbol: str) -> UnitDefinition:
    if symbol not in UNITS:
        raise UnknownUnit(f"[pouring] unknown unit '{symbol}'")
    return UNITS[symbol]


def dimension_of(symbol: str) -> str:
    return definition(symbol).dimension


@dataclass(frozen=True)
class Quantity:
    """
    A measurement, not a number.

    "about 5 mL", "5.00 mL", and "at least 5 mL" are different claims, and
    flattening them loses information a chemist deliberately recorded.
    """

    value: Decimal
    unit: str
    uncertainty: Optional[Decimal] = None
    qualifier: str = "exact"  # exact | approximate | lessThan | greaterThan

    def __post_init__(self) -> None:
        definition(self.unit)  # raises UnknownUnit

    @property
    def dimension(self) -> str:
        return dimension_of(self.unit)

    def to(self, unit: str) -> "Quantity":
        return convert(self, unit)

    def __str__(self) -> str:
        prefix = {
            "approximate": "~",
            "lessThan": "<",
            "greaterThan": ">",
            "exact": "",
        }[self.qualifier]
        tail = f" ± {self.uncertainty}" if self.uncertainty is not None else ""
        return f"{prefix}{self.value} {self.unit}{tail}"


def quantity(
    value: str | int | Decimal,
    unit: str,
    uncertainty: str | Decimal | None = None,
    qualifier: str = "exact",
) -> Quantity:
    """Builds a Quantity, taking value as a string to preserve precision."""
    return Quantity(
        value=Decimal(str(value)),
        unit=unit,
        uncertainty=None if uncertainty is None else Decimal(str(uncertainty)),
        qualifier=qualifier,
    )


def convert(source: Quantity, unit: str) -> Quantity:
    """
    Converts between units of the same dimension.

    Affine by construction, so Celsius and gauge pressure work without a
    special case at the call site.
    """
    target = definition(unit)
    origin = definition(source.unit)

    if origin.dimension != target.dimension:
        raise DimensionMismatch(
            f"[pouring] cannot convert {origin.dimension} "
            f"('{source.unit}') to {target.dimension} ('{unit}')"
        )

    base = source.value * origin.scale + origin.offset
    converted = (base - target.offset) / target.scale

    scaled_uncertainty = (
        None
        if source.uncertainty is None
        # Uncertainty is an interval width: it scales but never shifts.
        else source.uncertainty * origin.scale / target.scale
    )

    return Quantity(
        value=converted,
        unit=unit,
        uncertainty=scaled_uncertainty,
        qualifier=source.qualifier,
    )


def same_dimension(first: Quantity, second: Quantity) -> bool:
    return first.dimension == second.dimension


def add(first: Quantity, second: Quantity) -> Quantity:
    """
    Adds two quantities, in the units of the first.

    Rejects mismatched dimensions. Also rejects affine units, where addition
    is meaningless: 20 degC plus 20 degC is not 40 degC, because Celsius has
    an arbitrary zero. Convert to Kelvin first if a sum is genuinely intended.
    """
    if not same_dimension(first, second):
        raise DimensionMismatch(
            f"[pouring] cannot add {first.dimension} to {second.dimension}"
        )

    if definition(first.unit).offset != 0:
        raise DimensionMismatch(
            f"[pouring] '{first.unit}' has an arbitrary zero; addition is "
            f"not meaningful. Convert to an absolute unit first."
        )

    return Quantity(
        value=first.value + convert(second, first.unit).value,
        unit=first.unit,
    )
