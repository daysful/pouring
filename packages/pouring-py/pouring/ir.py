"""
The `pouring` IR, as data.

Mirrors about/notes/ir-specification.md. Everything here is plain, serialisable
structure — builders produce it, validators check it, simulators consume it. No
chemistry logic lives on these types.

Two target profiles share these entities and differ in how they are read:

    pouring:synthesis-v0   structural species, balance required, routes acyclic
    pouring:crn-v0         formal species allowed, balance N/A, cycles expected
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
from fractions import Fraction
from typing import Optional, Union

from .units import Quantity

SCHEMA_VERSION = "0.1.0"

SYNTHESIS = "pouring:synthesis-v0"
CRN = "pouring:crn-v0"

CHEMISTRY_V0 = "pouring:organic-covalent-v0"
NORMALIZE_V0 = "pouring:normalize-v0"

DETERMINISTIC = "deterministic"
STOCHASTIC = "stochastic"


# --- structure layer (spec section 4) --------------------------------------


@dataclass(frozen=True)
class Atom:
    id: str
    element: str
    charge: int = 0
    # Unpaired electrons. Spec 4.5: an atom short of its valence is an error,
    # never a silent radical — unpaired electrons must be declared.
    radical: int = 0
    isotope: Optional[int] = None
    label: Optional[str] = None

    def to_json(self) -> dict:
        out = {"id": self.id, "element": self.element, "charge": self.charge}
        if self.radical:
            out["radical"] = self.radical
        if self.isotope is not None:
            out["isotope"] = self.isotope
        return out


@dataclass(frozen=True)
class Bond:
    id: str
    source: str
    target: str
    order: int = 1

    def key(self) -> frozenset[str]:
        """Bonds are undirected: {a,b} and {b,a} are the same bond."""
        return frozenset((self.source, self.target))

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "from": self.source,
            "to": self.target,
            "order": self.order,
        }


@dataclass(frozen=True)
class TetrahedralCenter:
    """
    Spec 4.7. An ordered neighbour list plus a parity bit.

    Reading down the bond from `neighbors[0]` toward the centre, the remaining
    three neighbours appear either clockwise (+1) or anticlockwise (-1). This
    is what separates D-glucose from L-glucose, which are otherwise identical
    in every atom and bond.
    """

    atom: str
    neighbors: tuple[str, ...]
    parity: int

    def to_json(self) -> dict:
        return {
            "atom": self.atom,
            "neighbors": list(self.neighbors),
            "parity": self.parity,
        }


@dataclass(frozen=True)
class DoubleBondStereo:
    """Spec 4.7. E/Z, relative to one named neighbour on each end."""

    bond: str
    reference_source: str
    reference_target: str
    config: str  # "cis" | "trans"

    def to_json(self) -> dict:
        return {
            "bond": self.bond,
            "referenceFrom": self.reference_source,
            "referenceTo": self.reference_target,
            "config": self.config,
        }


@dataclass
class Molecule:
    """A structural species: a connected molecular graph."""

    id: str
    name: str
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    stereo_centers: list[TetrahedralCenter] = field(default_factory=list)
    stereo_bonds: list[DoubleBondStereo] = field(default_factory=list)

    kind = "molecule"

    def atom(self, atom_id: str) -> Atom:
        for candidate in self.atoms:
            if candidate.id == atom_id:
                return candidate
        raise KeyError(f"[pouring] no atom '{atom_id}' in '{self.name}'")

    def formula_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for atom in self.atoms:
            counts[atom.element] = counts.get(atom.element, 0) + 1
        return counts

    def total_charge(self) -> int:
        return sum(atom.charge for atom in self.atoms)

    def neighbors(self, atom_id: str) -> list[tuple[str, int]]:
        """Adjacent atom ids with the order of the connecting bond."""
        found = []
        for bond in self.bonds:
            if bond.source == atom_id:
                found.append((bond.target, bond.order))
            elif bond.target == atom_id:
                found.append((bond.source, bond.order))
        return found

    def to_json(self) -> dict:
        out = {
            "kind": self.kind,
            "id": self.id,
            "name": self.name,
            "atoms": [a.to_json() for a in self.atoms],
            "bonds": [b.to_json() for b in self.bonds],
        }
        if self.stereo_centers:
            out["stereoCenters"] = [s.to_json() for s in self.stereo_centers]
        if self.stereo_bonds:
            out["stereoBonds"] = [s.to_json() for s in self.stereo_bonds]
        return out


@dataclass
class FormalSpecies:
    """
    Identity without structure (spec 4.8).

    In a reaction network this is the normal case: the network is a program and
    its species are variables. They acquire structure only when the network is
    compiled to a physical realisation, at which point `binding` records what
    realised them.
    """

    id: str
    name: str
    binding: Optional[str] = None

    kind = "formalSpecies"

    def to_json(self) -> dict:
        out = {"kind": self.kind, "id": self.id, "name": self.name}
        if self.binding is not None:
            out["binding"] = self.binding
        return out


@dataclass(frozen=True)
class Component:
    entity: str
    ratio: Fraction = Fraction(1)

    def to_json(self) -> dict:
        return {
            "entity": self.entity,
            "ratio": {
                "numerator": str(self.ratio.numerator),
                "denominator": str(self.ratio.denominator),
            },
        }


@dataclass
class MultiComponentSpecies:
    """
    Spec 4.6. Salts, hydrates, solvates, co-crystals.

    Sodium chloride is sodium and chloride in 1:1 ratio, not a covalent
    molecule with a strange bond. Modelling it as one disconnected graph was
    the alternative, and it is wrong.
    """

    id: str
    name: str
    components: list[Component] = field(default_factory=list)

    kind = "multiComponentSpecies"

    def to_json(self) -> dict:
        return {
            "kind": self.kind,
            "id": self.id,
            "name": self.name,
            "components": [c.to_json() for c in self.components],
        }


Species = Union[Molecule, FormalSpecies, MultiComponentSpecies]


def is_structural(species: Species) -> bool:
    """True when atom-level checks (balance, formula, mass) can run."""
    return isinstance(species, (Molecule, MultiComponentSpecies))


# --- transformation layer (spec section 6) ---------------------------------


@dataclass(frozen=True)
class Participant:
    species: str
    coefficient: Fraction = Fraction(1)

    def to_json(self) -> dict:
        return {
            "species": self.species,
            "coefficient": {
                "numerator": str(self.coefficient.numerator),
                "denominator": str(self.coefficient.denominator),
            },
        }


@dataclass(frozen=True)
class Kinetics:
    """
    Spec 6.8. Optional under the synthesis target; required under CRN, where
    the rate constants *are* the program.
    """

    rate_constant: float
    rate_law: str = "massAction"

    def to_json(self) -> dict:
        return {"rateLaw": self.rate_law, "rateConstant": self.rate_constant}


@dataclass
class Conditions:
    """Spec 6.1. Typed quantities, not bare numbers."""

    temperature: Optional[Quantity] = None
    time: Optional[Quantity] = None
    pressure: Optional[Quantity] = None
    solvent: Optional[str] = None
    catalyst: Optional[str] = None

    def to_json(self) -> dict:
        out: dict = {}
        for name in ("temperature", "time", "pressure"):
            value = getattr(self, name)
            if value is not None:
                out[name] = {"value": str(value.value), "unit": value.unit}
        for name in ("solvent", "catalyst"):
            value = getattr(self, name)
            if value is not None:
                out[name] = value
        return out


@dataclass(frozen=True)
class Yield:
    """
    Spec 6.5. A single number is insufficient — a yield needs to say what it
    is a yield *of*, against what basis, and how it was determined.
    """

    product: str
    basis: str
    value: Decimal
    kind: str = "isolated"  # isolated | assay | calculated
    method: Optional[str] = None

    def to_json(self) -> dict:
        out = {
            "product": self.product,
            "basis": self.basis,
            "value": str(self.value),
            "kind": self.kind,
        }
        if self.method:
            out["method"] = self.method
        return out


# Spec 6.4 — not every recorded equation is meant to balance.
AUTHOR_SPECIFIED = "authorSpecified"
BALANCED = "balanced"
BALANCE_SUGGESTED = "balanceSuggested"
PARTIAL = "partial"
OBSERVATIONAL = "observational"


@dataclass
class Reaction:
    id: str
    name: str
    reactants: list[Participant]
    products: list[Participant]
    kinetics: Optional[Kinetics] = None
    conditions: Optional[Conditions] = None
    equation_status: str = BALANCED
    reaction_yield: Optional[Yield] = None
    # Spec 6.2: reactant atom id -> product atom id. Without it this is an
    # equation, not a transform, and cannot support retrosynthesis.
    atom_map: Optional[dict[str, str]] = None

    kind = "reactionEquation"

    def net_change(self) -> dict[str, Fraction]:
        """Stoichiometric change per species: products minus reactants."""
        delta: dict[str, Fraction] = {}
        for p in self.reactants:
            delta[p.species] = delta.get(p.species, Fraction(0)) - p.coefficient
        for p in self.products:
            delta[p.species] = delta.get(p.species, Fraction(0)) + p.coefficient
        return {s: d for s, d in delta.items() if d != 0}

    def to_json(self) -> dict:
        out = {
            "kind": self.kind,
            "id": self.id,
            "name": self.name,
            "reactants": [p.to_json() for p in self.reactants],
            "products": [p.to_json() for p in self.products],
            "equationStatus": self.equation_status,
        }
        if self.kinetics is not None:
            out["kinetics"] = self.kinetics.to_json()
        if self.conditions is not None:
            out["conditions"] = self.conditions.to_json()
        if self.reaction_yield is not None:
            out["yield"] = self.reaction_yield.to_json()
        if self.atom_map is not None:
            out["atomMap"] = dict(self.atom_map)
        return out


@dataclass
class ReactionNetwork:
    """
    Spec 6.7 — the CRN reading of a reaction set.

    No target and no termination requirement. Evaluated by setting initial
    amounts on `inputs` and observing `outputs` over time.
    """

    id: str
    name: str
    species: list[Species]
    reactions: list[Reaction]
    inputs: dict[str, float]
    outputs: list[str]
    semantics: Optional[str] = None

    kind = "reactionNetwork"

    def species_ids(self) -> list[str]:
        return [s.id for s in self.species]

    def to_json(self) -> dict:
        return {
            "kind": self.kind,
            "id": self.id,
            "name": self.name,
            "species": [s.to_json() for s in self.species],
            "reactions": [r.to_json() for r in self.reactions],
            "inputs": self.inputs,
            "outputs": self.outputs,
            "semantics": self.semantics,
        }


@dataclass
class Route:
    """
    Spec 6.6 — the synthesis reading of a reaction set. A plan: it has a
    target, it terminates, and it is acyclic.
    """

    id: str
    name: str
    target: str
    starting_materials: list[str]
    species: list[Species]
    reactions: list[Reaction]

    kind = "route"

    def to_json(self) -> dict:
        return {
            "kind": self.kind,
            "id": self.id,
            "name": self.name,
            "target": self.target,
            "startingMaterials": self.starting_materials,
            "species": [s.to_json() for s in self.species],
            "reactions": [r.to_json() for r in self.reactions],
        }


# --- document envelope (spec section 3) ------------------------------------


@dataclass
class Document:
    target_profile: str
    body: Union[ReactionNetwork, Route]
    chemistry_profile: str = CHEMISTRY_V0
    normalization_profile: str = NORMALIZE_V0
    schema_version: str = SCHEMA_VERSION

    def to_json(self) -> dict:
        return {
            "schemaVersion": self.schema_version,
            "targetProfile": self.target_profile,
            "chemistryProfile": self.chemistry_profile,
            "normalizationProfile": self.normalization_profile,
            "body": self.body.to_json(),
        }


# --- diagnostics (spec section 9) ------------------------------------------


@dataclass(frozen=True)
class Diagnostic:
    severity: str
    code: str
    message: str
    path: Optional[str] = None

    def __str__(self) -> str:
        where = f" [{self.path}]" if self.path else ""
        return f"{self.severity:7} {self.code}{where}: {self.message}"


ERROR = "error"
WARNING = "warning"
INFO = "info"
