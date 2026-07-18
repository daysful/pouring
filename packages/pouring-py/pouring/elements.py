"""
Element data, and the seam where chemistry authority lives.

The implementation decision (about/notes/implementation-language.md) is that
`pouring` must not own claims about chemistry — RDKit encodes decades of edge
cases that a hand-written table cannot. This module honours that by separating
two very different kinds of data:

**Atomic masses are facts.** Standard atomic weights are published numbers, not
judgements. Hard-coding them is safe.

**Permitted valences are a model.** Which valences an element may carry depends
on charge, bonding context, aromaticity, and documented special cases, and the
specification concedes (4.5) that its own rule is "not specified precisely
enough for two independent implementations to agree".

So valence goes behind `ValenceOracle`. The built-in implementation is a lint
for the restricted v0 profile and is labelled as such; `RDKitValenceOracle` is
the intended authority. Callers depend on the protocol, never on the table.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

# Standard atomic weights, IUPAC 2021 abridged.
MASSES: dict[str, float] = {
    "H": 1.008,
    "He": 4.0026,
    "Li": 6.94,
    "Be": 9.0122,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Fe": 55.845,
    "Cu": 63.546,
    "Zn": 65.38,
    "Br": 79.904,
    "I": 126.90,
}

# Neutral valences under `pouring:organic-covalent-v0`.
_NEUTRAL_VALENCES: dict[str, list[int]] = {
    "H": [1],
    "He": [0],
    "Li": [1],
    "Be": [2],
    "B": [3],
    "C": [4],
    "N": [3, 5],
    "O": [2],
    "F": [1],
    "Ne": [0],
    "Na": [1],
    "Mg": [2],
    "Al": [3],
    "Si": [4],
    "P": [3, 5],
    "S": [2, 4, 6],
    "Cl": [1],
    "Ar": [0],
    "K": [1],
    "Ca": [2],
    "Br": [1],
    "I": [1],
}

# Elements carrying lone pairs, for which a positive charge frees an extra
# bonding position (ammonium N takes four bonds, hydronium O takes three).
# For electron-deficient elements the opposite holds, and both a carbocation
# and a carbanion take three.
_LONE_PAIR_BEARING = {"N", "O", "F", "P", "S", "Cl", "Br", "I"}


# Exact masses for the isotopes that actually get labelled in practice.
# Not a complete nuclide table — a full one belongs with RDKit. Unlisted
# isotopes fall back to the standard atomic weight, which is announced by
# `is_known_isotope` so a caller can tell an exact mass from an approximation.
ISOTOPE_MASSES: dict[tuple[str, int], float] = {
    ("H", 1): 1.007825,
    ("H", 2): 2.014102,
    ("H", 3): 3.016049,
    ("C", 12): 12.000000,
    ("C", 13): 13.003355,
    ("C", 14): 14.003242,
    ("N", 14): 14.003074,
    ("N", 15): 15.000109,
    ("O", 16): 15.994915,
    ("O", 17): 16.999132,
    ("O", 18): 17.999160,
    ("S", 32): 31.972071,
    ("S", 34): 33.967867,
    ("Cl", 35): 34.968853,
    ("Cl", 37): 36.965903,
    ("Br", 79): 78.918338,
    ("Br", 81): 80.916291,
}


def is_known(element: str) -> bool:
    return element in MASSES


def is_known_isotope(element: str, isotope: int) -> bool:
    return (element, isotope) in ISOTOPE_MASSES


def mass(element: str, isotope: int | None = None) -> float:
    """
    Mass in g/mol: the exact isotope mass when one is specified and known,
    otherwise the standard atomic weight.

    The distinction matters. D2O is 20.03 g/mol against water's 18.02, and a
    labelling study that silently used the standard weight would report the
    wrong thing throughout.
    """
    if element not in MASSES:
        raise KeyError(f"[pouring] no mass for element '{element}'")

    if isotope is not None:
        exact = ISOTOPE_MASSES.get((element, isotope))
        if exact is not None:
            return exact

    return MASSES[element]


@runtime_checkable
class ValenceOracle(Protocol):
    """What `pouring` is allowed to ask about valence."""

    name: str

    def supports(self, element: str) -> bool:
        ...

    def permitted(self, element: str, charge: int) -> list[int]:
        """Valences the element may carry at this formal charge."""
        ...


class ProfileValenceOracle:
    """
    The `pouring:organic-covalent-v0` lint.

    A useful check for small covalent organic molecules and nothing more. It
    does not model hypervalence beyond the listed cases, organometallics,
    multi-centre bonding, or the context dependence of real valence. Anything
    it cannot represent must surface as `chemistry.unsupported` rather than a
    plausible-looking wrong answer.
    """

    name = "pouring:organic-covalent-v0"

    def supports(self, element: str) -> bool:
        return element in _NEUTRAL_VALENCES

    def permitted(self, element: str, charge: int) -> list[int]:
        if element not in _NEUTRAL_VALENCES:
            raise KeyError(f"[pouring] no valence model for '{element}'")

        neutral = _NEUTRAL_VALENCES[element]

        if charge == 0:
            return list(neutral)

        if element in _LONE_PAIR_BEARING:
            # A lone pair becomes a bonding position when charge is positive,
            # and a bonding position becomes a lone pair when negative.
            shifted = [v + charge for v in neutral]
        else:
            # Electron-deficient: either sign costs a bonding position.
            shifted = [v - abs(charge) for v in neutral]

        return [v for v in shifted if v >= 0]


class RDKitValenceOracle:
    """
    The intended authority, once RDKit is a dependency.

    Deliberately not implemented: writing a second hand-rolled table here and
    calling it RDKit would be worse than the honest gap. Install with
    `pip install -e '.[chemistry]'` and this becomes the real thing.
    """

    name = "rdkit"

    def __init__(self) -> None:
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError as error:
            raise ImportError(
                "[pouring] RDKit is not installed; valence authority is "
                "unavailable. Install with: pip install -e '.[chemistry]'"
            ) from error

    def supports(self, element: str) -> bool:  # pragma: no cover
        from rdkit import Chem

        return Chem.GetPeriodicTable().GetAtomicNumber(element) > 0

    def permitted(self, element: str, charge: int) -> list[int]:  # pragma: no cover
        from rdkit import Chem

        table = Chem.GetPeriodicTable()
        number = table.GetAtomicNumber(element)
        return [v - charge for v in table.GetValenceList(number) if v >= 0]


DEFAULT_ORACLE: ValenceOracle = ProfileValenceOracle()


def default_oracle() -> ValenceOracle:
    """
    Prefers RDKit when present, falls back to the profile lint.

    The fallback is announced through the `name` attribute so diagnostics can
    say which authority produced them — a valence verdict is only as good as
    the model behind it, and callers deserve to know which one answered.
    """
    try:
        return RDKitValenceOracle()
    except ImportError:
        return ProfileValenceOracle()
