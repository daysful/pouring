"""
Stoichiometric balancing (spec 6.3).

Balance is a nullspace problem over the integer atom-and-charge matrix,
`A nu = 0`. Solving it is not the whole problem, and the cases below are what
separate a balancer from a one-line linear algebra call:

  * **no solution** — the equation cannot be balanced at all;
  * **nullspace of dimension > 1** — the equation is underdetermined, and the
    answer is a *choice* rather than a result;
  * **sign constraints** — a mathematical solution that puts a negative
    coefficient on a species is not a chemical one;
  * **minimal integers** — 2:1:2 and 4:2:4 are the same balance, and only one
    is canonical.

Exact rational arithmetic throughout. Floating point makes minimal-integer
normalisation unreliable, because a coefficient that should be 3 arrives as
2.9999999999999996 and the whole answer is scaled off it.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Callable, Optional

from .ir import (
    BALANCED,
    Molecule,
    MultiComponentSpecies,
    Reaction,
    Species,
)

CHARGE_ROW = "__charge__"

Resolver = Callable[[str], Species]


@dataclass
class BalanceResult:
    status: str  # balanced | unbalanceable | underdetermined
    coefficients: Optional[dict[str, int]] = None
    freedom: int = 0
    detail: str = ""
    # Nullspace basis when the equation is underdetermined, so a caller can
    # inspect the alternatives rather than take one on trust.
    basis: Optional[list[list[Fraction]]] = None

    @property
    def ok(self) -> bool:
        return self.status == "balanced"


# --- atom tallies ----------------------------------------------------------


def species_counts(species: Species, resolve: Resolver) -> dict[str, Fraction]:
    """
    Atom counts plus a charge pseudo-element, so charge balances alongside
    mass rather than in a separate pass.
    """
    counts: dict[str, Fraction] = {}

    if isinstance(species, Molecule):
        for element, count in species.formula_counts().items():
            counts[element] = counts.get(element, Fraction(0)) + count
        counts[CHARGE_ROW] = Fraction(species.total_charge())
        return counts

    if isinstance(species, MultiComponentSpecies):
        for component in species.components:
            inner = species_counts(resolve(component.entity), resolve)
            for element, count in inner.items():
                counts[element] = (
                    counts.get(element, Fraction(0)) + component.ratio * count
                )
        return counts

    raise TypeError(
        f"[pouring] '{species.name}' has no structure; atom counts are "
        f"undefined for formal species"
    )


def _tally(participants, resolve: Resolver) -> dict[str, Fraction]:
    total: dict[str, Fraction] = {}
    for participant in participants:
        counts = species_counts(resolve(participant.species), resolve)
        for element, count in counts.items():
            total[element] = (
                total.get(element, Fraction(0)) + participant.coefficient * count
            )
    return {e: c for e, c in total.items() if c != 0}


def is_balanced(reaction: Reaction, resolve: Resolver) -> bool:
    return _tally(reaction.reactants, resolve) == _tally(
        reaction.products, resolve
    )


def imbalance(reaction: Reaction, resolve: Resolver) -> tuple[dict, dict]:
    return _tally(reaction.reactants, resolve), _tally(reaction.products, resolve)


# --- exact linear algebra --------------------------------------------------


def _rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    """Reduced row echelon form over the rationals. Exact, no pivoting noise."""
    work = [row[:] for row in matrix]
    if not work:
        return work, []

    rows, columns = len(work), len(work[0])
    pivots: list[int] = []
    row = 0

    for column in range(columns):
        candidate = next(
            (r for r in range(row, rows) if work[r][column] != 0), None
        )
        if candidate is None:
            continue

        work[row], work[candidate] = work[candidate], work[row]

        divisor = work[row][column]
        work[row] = [value / divisor for value in work[row]]

        for other in range(rows):
            if other != row and work[other][column] != 0:
                factor = work[other][column]
                work[other] = [
                    a - factor * b for a, b in zip(work[other], work[row])
                ]

        pivots.append(column)
        row += 1
        if row == rows:
            break

    return work, pivots


def _nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    if not matrix:
        return []

    reduced, pivots = _rref(matrix)
    columns = len(matrix[0])
    free = [c for c in range(columns) if c not in pivots]

    basis: list[list[Fraction]] = []
    for column in free:
        vector = [Fraction(0)] * columns
        vector[column] = Fraction(1)
        for index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[index][column]
        basis.append(vector)

    return basis


def _minimal_integers(vector: list[Fraction]) -> Optional[list[int]]:
    """
    Scales a rational vector to the smallest all-positive integer vector.

    Returns None when the vector mixes signs — mathematically a solution,
    chemically not one, since a species cannot participate negatively.
    """
    nonzero = [value for value in vector if value != 0]
    if not nonzero:
        return None

    if all(value < 0 for value in nonzero):
        vector = [-value for value in vector]
        nonzero = [-value for value in nonzero]

    if any(value < 0 for value in nonzero):
        return None

    multiplier = 1
    for value in vector:
        multiplier = multiplier * value.denominator // gcd(
            multiplier, value.denominator
        )

    scaled = [int(value * multiplier) for value in vector]

    divisor = 0
    for value in scaled:
        divisor = gcd(divisor, abs(value))

    if divisor > 1:
        scaled = [value // divisor for value in scaled]

    return scaled


# --- the balancer ----------------------------------------------------------


def balance(reaction: Reaction, resolve: Resolver) -> BalanceResult:
    """
    Solves for stoichiometric coefficients.

    Reactants are entered positive and products negative, so a nullspace
    vector with all-positive entries is exactly a chemically valid balance.
    """
    species_order = [p.species for p in reaction.reactants] + [
        p.species for p in reaction.products
    ]
    signs = [1] * len(reaction.reactants) + [-1] * len(reaction.products)

    try:
        columns = [
            species_counts(resolve(species_id), resolve)
            for species_id in species_order
        ]
    except TypeError as error:
        return BalanceResult(status="unbalanceable", detail=str(error))

    elements = sorted({element for column in columns for element in column})

    matrix = [
        [
            sign * column.get(element, Fraction(0))
            for column, sign in zip(columns, signs)
        ]
        for element in elements
    ]

    basis = _nullspace(matrix)

    if not basis:
        return BalanceResult(
            status="unbalanceable",
            detail="no non-trivial solution; these species cannot balance",
        )

    if len(basis) > 1:
        # More than one independent solution: the equation does not determine
        # its own coefficients, and picking one silently would be a guess.
        #
        # This is commoner than it sounds. Aspirin's synthesis has a
        # two-dimensional nullspace, because its hydrogen row is exactly twice
        # its oxygen row. The 1:1:1:1 a chemist writes is correct, and so is
        # 0:10:2:11 — which consumes no salicylic acid and still produces
        # aspirin. Atom counting cannot tell them apart; only chemistry can.
        return BalanceResult(
            status="underdetermined",
            freedom=len(basis),
            basis=basis,
            detail=(
                f"nullspace has dimension {len(basis)}; the coefficients are "
                f"a choice, not a result"
            ),
        )

    integers = _minimal_integers(basis[0])

    if integers is None:
        return BalanceResult(
            status="unbalanceable",
            detail="only solutions require a negative coefficient",
        )

    coefficients: dict[str, int] = {}
    for species_id, value in zip(species_order, integers):
        # A species on both sides accumulates; net participation is what
        # balances.
        coefficients[species_id] = value

    return BalanceResult(status="balanced", coefficients=coefficients)


def balanced_equation(reaction: Reaction, resolve: Resolver) -> str:
    """Renders the solved equation, for demos and diagnostics."""
    result = balance(reaction, resolve)
    if not result.ok or result.coefficients is None:
        return f"<{result.status}: {result.detail}>"

    def side(participants) -> str:
        parts = []
        for participant in participants:
            count = result.coefficients[participant.species]
            name = resolve(participant.species).name
            parts.append(name if count == 1 else f"{count} {name}")
        return " + ".join(parts)

    return f"{side(reaction.reactants)} -> {side(reaction.products)}"


def suggested_coefficients(
    reaction: Reaction, resolve: Resolver
) -> Optional[dict[str, int]]:
    """Coefficients a `balanceSuggested` equation should adopt."""
    result = balance(reaction, resolve)
    return result.coefficients if result.ok else None


__all__ = [
    "BalanceResult",
    "CHARGE_ROW",
    "balance",
    "balanced_equation",
    "imbalance",
    "is_balanced",
    "species_counts",
    "suggested_coefficients",
]
