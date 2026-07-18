"""
Canonical labelling and content-addressed identity (spec 10.2).

Deterministic atom ids do not by themselves give chemical identity. This module
delivers the part that is honestly achievable now — a canonical form for the
*graph*, invariant under atom ordering, bond direction, and Kekulé choice — and
is explicit about the part that is not.

The algorithm is Morgan-style: refine atom invariants against neighbourhoods
until the partition stabilises, then break residual ties by trying each
candidate and keeping the lexicographically smallest encoding. For the
symmetric molecules where ties actually arise (benzene), every choice yields
the same encoding by symmetry, which is what makes the result canonical rather
than merely deterministic.

**Kekulé normalisation.** Benzene drawn 1-2-1-2-1-2 and drawn 2-1-2-1-2-1 are
the same molecule. Hashing bond orders directly would make them different, so
rings whose orders alternate are encoded with a symbolic aromatic bond, the
same device SMILES uses. The ring perception behind this is a narrow heuristic
for the v0 profile, not a general aromaticity model — spec 4.2 notes that
aromaticity is model-dependent and RDKit supports several.

**What this is not.** Chemical identity additionally requires policies for
resonance, tautomerism, protonation state, and salt disconnection. None are
implemented. A hash from this module means "the same annotated graph", which is
necessary for identity and not sufficient.
"""

from __future__ import annotations

import hashlib
from itertools import permutations

from .ir import CHEMISTRY_V0, NORMALIZE_V0, SCHEMA_VERSION, Molecule

MAX_RING_SIZE = 10
AROMATIC = "a"


# --- ring perception -------------------------------------------------------


def find_rings(molecule: Molecule) -> list[list[str]]:
    """
    Cycles up to MAX_RING_SIZE, as lists of atom ids.

    Sufficient for the small molecules the v0 profile targets. A full smallest
    set of smallest rings is a harder problem and is not attempted.
    """
    adjacency: dict[str, list[str]] = {atom.id: [] for atom in molecule.atoms}
    for bond in molecule.bonds:
        if bond.source in adjacency and bond.target in adjacency:
            adjacency[bond.source].append(bond.target)
            adjacency[bond.target].append(bond.source)

    seen: set[frozenset[str]] = set()
    rings: list[list[str]] = []

    def walk(start: str, current: str, path: list[str]) -> None:
        if len(path) > MAX_RING_SIZE:
            return
        for neighbor in adjacency[current]:
            if neighbor == start and len(path) >= 3:
                signature = frozenset(path)
                if signature not in seen:
                    seen.add(signature)
                    rings.append(list(path))
            elif neighbor not in path and neighbor > start:
                # `neighbor > start` keeps each cycle from being rediscovered
                # from every one of its members.
                walk(start, neighbor, path + [neighbor])

    for atom in molecule.atoms:
        walk(atom.id, atom.id, [atom.id])

    return rings


def _bond_between(molecule: Molecule, first: str, second: str):
    for bond in molecule.bonds:
        if bond.key() == frozenset((first, second)):
            return bond
    return None


def perceive_aromatic_bonds(molecule: Molecule) -> set[str]:
    """
    Bond ids belonging to a ring whose orders alternate single/double.

    A deliberately narrow heuristic: even-sized ring, orders strictly
    alternating between 1 and 2. It captures benzene and the Kekulé ambiguity
    that motivates this, and claims nothing more.
    """
    aromatic: set[str] = set()

    for ring in find_rings(molecule):
        size = len(ring)
        if size % 2 != 0:
            continue

        bonds = []
        for index, atom_id in enumerate(ring):
            bond = _bond_between(molecule, atom_id, ring[(index + 1) % size])
            if bond is None:
                break
            bonds.append(bond)

        if len(bonds) != size:
            continue

        orders = [bond.order for bond in bonds]
        alternating = all(
            {orders[i], orders[(i + 1) % size]} == {1, 2} for i in range(size)
        )

        if alternating:
            aromatic.update(bond.id for bond in bonds)

    return aromatic


# --- canonical labelling ---------------------------------------------------


def _initial_invariant(molecule: Molecule, atom_id: str, aromatic: set[str]) -> str:
    atom = molecule.atom(atom_id)
    degree = len(molecule.neighbors(atom_id))
    order_sum = sum(order for _, order in molecule.neighbors(atom_id))
    return "|".join(
        [
            atom.element,
            str(atom.charge),
            str(atom.radical),
            str(atom.isotope or 0),
            str(degree),
            str(order_sum),
        ]
    )


def _bond_symbol(molecule: Molecule, first: str, second: str, aromatic: set[str]) -> str:
    bond = _bond_between(molecule, first, second)
    if bond is None:
        return "?"
    return AROMATIC if bond.id in aromatic else str(bond.order)


def _refine(
    molecule: Molecule,
    invariants: dict[str, str],
    aromatic: set[str],
) -> dict[str, str]:
    """One refinement round: fold each atom's neighbourhood into its label."""
    updated: dict[str, str] = {}
    for atom in molecule.atoms:
        neighborhood = sorted(
            f"{_bond_symbol(molecule, atom.id, neighbor, aromatic)}"
            f"{invariants[neighbor]}"
            for neighbor, _ in molecule.neighbors(atom.id)
        )
        digest = hashlib.sha256(
            (invariants[atom.id] + "(" + ",".join(neighborhood) + ")").encode()
        ).hexdigest()[:16]
        updated[atom.id] = digest
    return updated


def _stabilise(
    molecule: Molecule,
    invariants: dict[str, str],
    aromatic: set[str],
) -> dict[str, str]:
    previous_classes = -1
    for _ in range(len(molecule.atoms) + 1):
        classes = len(set(invariants.values()))
        if classes == previous_classes:
            break
        previous_classes = classes
        invariants = _refine(molecule, invariants, aromatic)
    return invariants


def canonical_order(molecule: Molecule) -> list[str]:
    """
    Atom ids in canonical order.

    Refinement first; where it leaves atoms indistinguishable, each candidate
    is tried and the lexicographically smallest encoding wins.
    """
    aromatic = perceive_aromatic_bonds(molecule)

    invariants = {
        atom.id: _initial_invariant(molecule, atom.id, aromatic)
        for atom in molecule.atoms
    }
    invariants = _stabilise(molecule, invariants, aromatic)

    fixed: dict[str, str] = {}

    def resolve(current: dict[str, str]) -> list[str]:
        grouped: dict[str, list[str]] = {}
        for atom_id, label in current.items():
            grouped.setdefault(label, []).append(atom_id)

        tied = [ids for ids in grouped.values() if len(ids) > 1]

        if not tied:
            return sorted(current, key=lambda a: (current[a], a))

        # Break the smallest tie first: fewest branches for the same effect.
        smallest = min(tied, key=len)
        best: list[str] | None = None

        for candidate in sorted(smallest):
            attempt = dict(current)
            attempt[candidate] = "!" + attempt[candidate]
            attempt = _stabilise(molecule, attempt, aromatic)
            order = resolve(attempt)
            encoded = _encode_with_order(molecule, order, aromatic)
            if best is None or encoded < _encode_with_order(
                molecule, best, aromatic
            ):
                best = order

        return best if best is not None else sorted(current)

    return resolve(invariants)


# --- canonical encoding ----------------------------------------------------


def _permutation_sign(source: tuple[str, ...], target: tuple[str, ...]) -> int:
    """
    Sign of the permutation taking `source` to `target`, by inversion count.

    Stereo parity is declared against one neighbour ordering; expressing it in
    canonical order means composing with the permutation between them.
    """
    index = {value: position for position, value in enumerate(target)}
    sequence = [index[value] for value in source]

    inversions = sum(
        1
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
        if sequence[i] > sequence[j]
    )
    return 1 if inversions % 2 == 0 else -1


def _encode_with_order(
    molecule: Molecule,
    order: list[str],
    aromatic: set[str],
) -> str:
    position = {atom_id: index for index, atom_id in enumerate(order)}

    atoms = []
    for atom_id in order:
        atom = molecule.atom(atom_id)
        atoms.append(
            f"{atom.element},{atom.charge},{atom.radical},{atom.isotope or 0}"
        )

    bonds = []
    for bond in molecule.bonds:
        if bond.source not in position or bond.target not in position:
            continue
        low, high = sorted((position[bond.source], position[bond.target]))
        symbol = AROMATIC if bond.id in aromatic else str(bond.order)
        bonds.append(f"{low}-{high}:{symbol}")
    bonds.sort()

    stereo = []
    for center in molecule.stereo_centers:
        if center.atom not in position:
            continue
        canonical_neighbors = tuple(
            sorted(center.neighbors, key=lambda a: position.get(a, -1))
        )
        sign = _permutation_sign(tuple(center.neighbors), canonical_neighbors)
        stereo.append(f"{position[center.atom]}:{center.parity * sign}")
    stereo.sort()

    bond_stereo = []
    for entry in molecule.stereo_bonds:
        bond = next((b for b in molecule.bonds if b.id == entry.bond), None)
        if bond is None:
            continue
        low, high = sorted((position[bond.source], position[bond.target]))
        references = sorted(
            (
                position.get(entry.reference_source, -1),
                position.get(entry.reference_target, -1),
            )
        )
        bond_stereo.append(f"{low}-{high}:{references}:{entry.config}")
    bond_stereo.sort()

    return "/".join(
        [
            "A:" + ";".join(atoms),
            "B:" + ";".join(bonds),
            "S:" + ";".join(stereo),
            "D:" + ";".join(bond_stereo),
        ]
    )


def canonical_encoding(molecule: Molecule) -> str:
    """Ordering-independent string form of the annotated graph."""
    aromatic = perceive_aromatic_bonds(molecule)
    return _encode_with_order(molecule, canonical_order(molecule), aromatic)


def content_hash(
    molecule: Molecule,
    chemistry_profile: str = CHEMISTRY_V0,
    normalization_profile: str = NORMALIZE_V0,
    schema_version: str = SCHEMA_VERSION,
) -> str:
    """
    Spec 10.2:

        hash(schemaVersion, chemistryProfile, normalizationProfile,
             canonicalStructureEncoding)

    The profiles are inside the hash on purpose. Hashing the graph alone would
    make two documents equal that assert different things, because the same
    graph under a different chemistry model is not the same claim.
    """
    payload = "\n".join(
        [
            schema_version,
            chemistry_profile,
            normalization_profile,
            canonical_encoding(molecule),
        ]
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def same_structure(first: Molecule, second: Molecule) -> bool:
    return canonical_encoding(first) == canonical_encoding(second)
