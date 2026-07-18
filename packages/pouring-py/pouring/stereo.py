"""
Stereochemistry (spec 4.7).

The specification calls this the largest gap and disqualifying for exact
chemical identity, because D-glucose and L-glucose have identical atoms, bonds,
and formula and differ only here. This module closes the first two items of
that section — tetrahedral parity and E/Z configuration — and detects centres
that are potentially chiral but unannotated.

**Substituent distinctness is a heuristic here, not CIP.** Deciding whether a
carbon's four substituents differ is done by comparing a bounded fingerprint of
each branch. Full Cahn-Ingold-Prelog priority is a substantial ruleset with
documented edge cases, and belongs behind the RDKit reference rather than being
approximated here. The fingerprint is conservative in the direction that
matters: it may miss a stereocentre, and reports what it did not decide.

**Not implemented:** axial chirality (allenes, atropisomerism), square-planar
and octahedral coordination, and enhanced stereo groups. Spec 4.7 places these
outside the v0 profile, so they raise `chemistry.unsupported` rather than being
silently dropped.
"""

from __future__ import annotations

from .ir import INFO, WARNING, Diagnostic, Molecule, TetrahedralCenter

MAX_BRANCH_DEPTH = 6


def _branch_signature(molecule: Molecule, start: str, blocked: str) -> str:
    """
    A bounded fingerprint of everything reachable from `start` without passing
    back through `blocked`.

    Two substituents with equal signatures are treated as indistinguishable.
    Bounded depth means deep-but-different branches can collide, which is why
    this is stated as conservative rather than exact.
    """
    layers: list[str] = []
    frontier = {start}
    visited = {start, blocked}

    for _ in range(MAX_BRANCH_DEPTH):
        if not frontier:
            break

        described = sorted(
            f"{molecule.atom(atom_id).element}"
            f"{molecule.atom(atom_id).charge:+d}"
            for atom_id in frontier
        )
        layers.append(",".join(described))

        following: set[str] = set()
        for atom_id in frontier:
            for neighbor, order in molecule.neighbors(atom_id):
                if neighbor not in visited:
                    visited.add(neighbor)
                    following.add(neighbor)
        frontier = following

    return "/".join(layers)


def potential_stereocenters(molecule: Molecule) -> list[str]:
    """
    Atoms with four mutually distinct substituents.

    Restricted to tetrahedral carbon-like centres: exactly four single-bonded
    neighbours. Centres with a double bond are not tetrahedral and are handled
    by E/Z instead.
    """
    candidates: list[str] = []

    for atom in molecule.atoms:
        neighbors = molecule.neighbors(atom.id)

        if len(neighbors) != 4:
            continue
        if any(order != 1 for _, order in neighbors):
            continue

        signatures = {
            _branch_signature(molecule, neighbor, atom.id)
            for neighbor, _ in neighbors
        }

        if len(signatures) == 4:
            candidates.append(atom.id)

    return candidates


def annotated_centers(molecule: Molecule) -> set[str]:
    return {center.atom for center in molecule.stereo_centers}


def check_stereo(
    molecule: Molecule,
    severity: str = WARNING,
) -> list[Diagnostic]:
    """
    Reports potentially chiral centres carrying no annotation.

    Severity is a parameter because spec 4.7 makes it contextual: for an exact
    active ingredient an unannotated centre is an error, for a declared
    racemate it is intentional, for an intermediate where stereochemistry is
    irrelevant it is acceptable, and for an imported historical record it is
    simply unknown. This module cannot know which, so the caller decides.
    """
    annotated = annotated_centers(molecule)

    return [
        Diagnostic(
            severity,
            "stereo.unspecified",
            f"'{atom_id}' has four distinct substituents but no configuration; "
            f"this structure does not distinguish its enantiomers",
            path=atom_id,
        )
        for atom_id in potential_stereocenters(molecule)
        if atom_id not in annotated
    ]


def invert(center: TetrahedralCenter) -> TetrahedralCenter:
    """The opposite configuration at one centre."""
    return TetrahedralCenter(
        atom=center.atom,
        neighbors=center.neighbors,
        parity=-center.parity,
    )


def mirror(molecule: Molecule, name: str | None = None) -> Molecule:
    """
    The enantiomer: every tetrahedral centre inverted.

    Used to prove the representation actually distinguishes them — a mirrored
    molecule must have the same formula and mass and a different canonical
    form, and if it does not, the stereo layer is decorative.
    """
    return Molecule(
        id=f"{molecule.id}-ent",
        name=name or f"ent-{molecule.name}",
        atoms=list(molecule.atoms),
        bonds=list(molecule.bonds),
        stereo_centers=[invert(center) for center in molecule.stereo_centers],
        stereo_bonds=list(molecule.stereo_bonds),
    )


def stereo_summary(molecule: Molecule) -> str:
    potential = potential_stereocenters(molecule)
    annotated = annotated_centers(molecule)
    unspecified = [atom_id for atom_id in potential if atom_id not in annotated]

    return (
        f"{len(potential)} potential centre(s), "
        f"{len(annotated & set(potential))} annotated, "
        f"{len(unspecified)} unspecified"
    )
