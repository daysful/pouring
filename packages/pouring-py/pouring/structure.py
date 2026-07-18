"""
Structure-layer analysis and checking (spec section 4).

Everything here operates on the molecular graph: formula, mass, connectivity,
rings, and the valence rule of 4.5. Valence verdicts come from a
`ValenceOracle` rather than from this module, because which valences an element
may carry is a model and not a fact — see elements.py.
"""

from __future__ import annotations

from decimal import Decimal

from .elements import ValenceOracle, is_known, mass
from .ir import (
    ERROR,
    Diagnostic,
    Molecule,
    MultiComponentSpecies,
    Species,
    WARNING,
)


# --- formula and mass ------------------------------------------------------


def hill_formula(counts: dict[str, int]) -> str:
    """
    Hill notation: carbon first, hydrogen second, everything else
    alphabetical. When there is no carbon, everything is alphabetical.
    """
    if not counts:
        return ""

    remaining = dict(counts)
    ordered: list[tuple[str, int]] = []

    if "C" in remaining:
        ordered.append(("C", remaining.pop("C")))
        if "H" in remaining:
            ordered.append(("H", remaining.pop("H")))

    ordered.extend(sorted(remaining.items()))

    return "".join(
        element if count == 1 else f"{element}{count}" for element, count in ordered
    )


def formula(molecule: Molecule) -> str:
    return hill_formula(molecule.formula_counts())


def molar_mass(molecule: Molecule) -> Decimal:
    """
    Mass in g/mol, honouring declared isotopes.

    Labelled atoms use their exact isotope mass where known; everything else
    uses the standard atomic weight.
    """
    total = Decimal(0)
    for atom in molecule.atoms:
        total += Decimal(str(mass(atom.element, atom.isotope)))
    return total


def component_counts(
    species: MultiComponentSpecies,
    resolve,
) -> dict[str, Decimal]:
    """
    Atom counts for a multi-component species, weighted by component ratio.

    Ratios may be fractional (a hemihydrate is 1:0.5), so counts are Decimal
    rather than int.
    """
    counts: dict[str, Decimal] = {}
    for component in species.components:
        entity = resolve(component.entity)
        ratio = Decimal(component.ratio.numerator) / Decimal(
            component.ratio.denominator
        )
        for element, count in entity.formula_counts().items():
            counts[element] = counts.get(element, Decimal(0)) + ratio * count
    return counts


# --- connectivity ----------------------------------------------------------


def connected_components(molecule: Molecule) -> list[set[str]]:
    adjacency: dict[str, set[str]] = {atom.id: set() for atom in molecule.atoms}
    for bond in molecule.bonds:
        if bond.source in adjacency and bond.target in adjacency:
            adjacency[bond.source].add(bond.target)
            adjacency[bond.target].add(bond.source)

    unvisited = set(adjacency)
    components: list[set[str]] = []

    while unvisited:
        start = unvisited.pop()
        component = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if neighbor not in component:
                    component.add(neighbor)
                    unvisited.discard(neighbor)
                    stack.append(neighbor)
        components.append(component)

    return components


def ring_count(molecule: Molecule) -> int:
    """
    Cyclomatic number: bonds - atoms + components.

    Benzene gives 1, which is the point of storing structure as a graph rather
    than a tree — a tree cannot have a ring at all.
    """
    return (
        len(molecule.bonds)
        - len(molecule.atoms)
        + len(connected_components(molecule))
    )


# --- valence ---------------------------------------------------------------


def used_valence(molecule: Molecule, atom_id: str) -> int:
    """
    Spec 4.5:  used(a) = sum of incident bond orders + radical electrons.

    Ports are absent because templates are not implemented; a Molecule is
    always fully elaborated.
    """
    atom = molecule.atom(atom_id)
    bonded = sum(order for _, order in molecule.neighbors(atom_id))
    return bonded + atom.radical


def check_structure(
    molecule: Molecule,
    oracle: ValenceOracle,
) -> list[Diagnostic]:
    """
    All structure-layer checks for one molecule.

    Ordering matters: graph integrity is established before valence, because
    valence computed over dangling bonds would produce confident nonsense.
    """
    found: list[Diagnostic] = []

    if not molecule.atoms:
        return [
            Diagnostic(
                ERROR, "molecule.empty", f"'{molecule.name}' has no atoms",
                path=molecule.id,
            )
        ]

    known_atoms = {atom.id for atom in molecule.atoms}
    graph_is_sound = True

    for atom in molecule.atoms:
        if not is_known(atom.element):
            found.append(
                Diagnostic(
                    ERROR,
                    "element.unknown",
                    f"'{atom.element}' is not a known element",
                    path=atom.id,
                )
            )
            graph_is_sound = False
        elif not oracle.supports(atom.element):
            found.append(
                Diagnostic(
                    ERROR,
                    "chemistry.unsupported",
                    f"'{atom.element}' is outside the {oracle.name} profile",
                    path=atom.id,
                )
            )
            graph_is_sound = False

    seen_pairs: set[frozenset[str]] = set()
    for bond in molecule.bonds:
        if bond.source not in known_atoms or bond.target not in known_atoms:
            found.append(
                Diagnostic(
                    ERROR,
                    "bond.dangling",
                    f"bond '{bond.id}' references an unknown atom",
                    path=bond.id,
                )
            )
            graph_is_sound = False
            continue

        if bond.source == bond.target:
            found.append(
                Diagnostic(
                    ERROR,
                    "bond.selfLoop",
                    f"bond '{bond.id}' joins '{bond.source}' to itself",
                    path=bond.id,
                )
            )
            graph_is_sound = False
            continue

        if bond.key() in seen_pairs:
            found.append(
                Diagnostic(
                    ERROR,
                    "bond.duplicate",
                    f"bond '{bond.id}' repeats a pair already bonded; raise "
                    f"the bond order instead",
                    path=bond.id,
                )
            )
            graph_is_sound = False
        seen_pairs.add(bond.key())

        if bond.order not in (1, 2, 3):
            found.append(
                Diagnostic(
                    ERROR,
                    "chemistry.unsupported",
                    f"bond order {bond.order} is outside the profile",
                    path=bond.id,
                )
            )
            graph_is_sound = False

    components = connected_components(molecule)
    if len(components) > 1:
        found.append(
            Diagnostic(
                ERROR,
                "molecule.disconnected",
                f"'{molecule.name}' has {len(components)} disconnected parts; "
                f"use a multi-component species instead",
                path=molecule.id,
            )
        )

    if not graph_is_sound:
        return found

    for atom in molecule.atoms:
        used = used_valence(molecule, atom.id)
        permitted = oracle.permitted(atom.element, atom.charge)

        if not permitted:
            continue

        if used in permitted:
            continue

        if used > max(permitted):
            found.append(
                Diagnostic(
                    ERROR,
                    "valence.exceeded",
                    f"{atom.element} '{atom.id}' carries {used} bond-order "
                    f"units, above the maximum {max(permitted)} "
                    f"({oracle.name})",
                    path=atom.id,
                )
            )
        else:
            short = min(v for v in permitted if v > used) - used
            found.append(
                Diagnostic(
                    ERROR,
                    "valence.floating",
                    f"{atom.element} '{atom.id}' carries {used} of "
                    f"{permitted}; {short} unsatisfied. Declare a radical if "
                    f"this is intended ({oracle.name})",
                    path=atom.id,
                )
            )

    found.extend(_check_stereo_wellformed(molecule))

    return found


def _check_stereo_wellformed(molecule: Molecule) -> list[Diagnostic]:
    """
    Stereo annotations must refer to real atoms and real bonds. Whether they
    are chemically *complete* is a separate question — see stereo.py.
    """
    found = []
    known_atoms = {atom.id for atom in molecule.atoms}
    known_bonds = {bond.id: bond for bond in molecule.bonds}

    for center in molecule.stereo_centers:
        if center.atom not in known_atoms:
            found.append(
                Diagnostic(
                    ERROR,
                    "stereo.dangling",
                    f"stereocentre references unknown atom '{center.atom}'",
                    path=center.atom,
                )
            )
            continue

        actual = {neighbor for neighbor, _ in molecule.neighbors(center.atom)}
        declared = set(center.neighbors)

        if declared != actual:
            found.append(
                Diagnostic(
                    ERROR,
                    "stereo.neighborMismatch",
                    f"stereocentre '{center.atom}' lists {sorted(declared)} "
                    f"but is bonded to {sorted(actual)}",
                    path=center.atom,
                )
            )

        if center.parity not in (1, -1):
            found.append(
                Diagnostic(
                    ERROR,
                    "stereo.invalidParity",
                    f"parity must be +1 or -1, got {center.parity}",
                    path=center.atom,
                )
            )

    for stereo_bond in molecule.stereo_bonds:
        bond = known_bonds.get(stereo_bond.bond)
        if bond is None:
            found.append(
                Diagnostic(
                    ERROR,
                    "stereo.dangling",
                    f"stereo bond references unknown bond '{stereo_bond.bond}'",
                    path=stereo_bond.bond,
                )
            )
            continue

        if bond.order != 2:
            found.append(
                Diagnostic(
                    WARNING,
                    "stereo.notDoubleBond",
                    f"E/Z configuration declared on bond '{bond.id}' of order "
                    f"{bond.order}",
                    path=bond.id,
                )
            )

        if stereo_bond.config not in ("cis", "trans"):
            found.append(
                Diagnostic(
                    ERROR,
                    "stereo.invalidConfig",
                    f"configuration must be cis or trans, got "
                    f"'{stereo_bond.config}'",
                    path=stereo_bond.bond,
                )
            )

    return found


def describe(species: Species) -> str:
    """A short human-readable summary, for demos and diagnostics."""
    if isinstance(species, Molecule):
        return (
            f"{species.name}: {formula(species)}  "
            f"{molar_mass(species):.3f} g/mol  "
            f"{ring_count(species)} ring(s)"
        )
    return f"{species.name}: (no structure)"
