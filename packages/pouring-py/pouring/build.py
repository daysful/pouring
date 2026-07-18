"""
Construction of molecular graphs.

Atoms are instances, bonds are nets, and because bonds are declared
independently of the order atoms appear, rings close naturally. This is the
netlist idea from the specification made concrete: `ring()` is the operation a
nested expression tree cannot express.

`fill()` is the elaboration step. Chemists write implicit hydrogens; the IR
stores every atom (spec 1.3), so the builder saturates open valences on the way
out. What it saturates *to* comes from the valence oracle, not from this
module.
"""

from __future__ import annotations

from typing import Callable, Optional

from .elements import ValenceOracle, default_oracle
from .ir import Atom, Bond, DoubleBondStereo, Molecule, TetrahedralCenter


class MoleculeBuilder:
    def __init__(self, name: str, oracle: Optional[ValenceOracle] = None) -> None:
        self.name = name
        self._oracle = oracle or default_oracle()
        self._atoms: list[Atom] = []
        self._bonds: list[Bond] = []
        self._centers: list[TetrahedralCenter] = []
        self._stereo_bonds: list[DoubleBondStereo] = []
        self._counter = 0

    def _next(self, prefix: str) -> str:
        self._counter += 1
        return f"{prefix}{self._counter}"

    # --- atoms and bonds ---------------------------------------------------

    def atom(
        self,
        element: str,
        charge: int = 0,
        radical: int = 0,
        isotope: Optional[int] = None,
        label: Optional[str] = None,
    ) -> str:
        atom = Atom(
            id=label or self._next("a"),
            element=element,
            charge=charge,
            radical=radical,
            isotope=isotope,
        )
        self._atoms.append(atom)
        return atom.id

    def atoms(self, element: str, count: int) -> list[str]:
        return [self.atom(element) for _ in range(count)]

    def bond(self, source: str, target: str, order: int = 1) -> str:
        if source == target:
            raise ValueError(
                f"[pouring] '{self.name}': cannot bond '{source}' to itself"
            )
        bond = Bond(id=self._next("b"), source=source, target=target, order=order)
        self._bonds.append(bond)
        return bond.id

    def chain(self, atoms: list[str], orders: Optional[list[int]] = None) -> list[str]:
        """An open sequence: n atoms, n-1 bonds."""
        return [
            self.bond(atom, atoms[index + 1], orders[index] if orders else 1)
            for index, atom in enumerate(atoms[:-1])
        ]

    def ring(self, atoms: list[str], orders: Optional[list[int]] = None) -> list[str]:
        """A closed cycle: n atoms, n bonds."""
        if len(atoms) < 3:
            raise ValueError(f"[pouring] '{self.name}': a ring needs 3+ atoms")

        size = len(atoms)
        return [
            self.bond(
                atom,
                atoms[(index + 1) % size],
                orders[index] if orders else 1,
            )
            for index, atom in enumerate(atoms)
        ]

    def group(self, anchor: str, element: str, order: int = 1) -> str:
        """Attaches one new atom to an existing one, returning the new atom."""
        created = self.atom(element)
        self.bond(anchor, created, order)
        return created

    def hydroxyl(self, anchor: str) -> str:
        """-OH. The hydrogen arrives via fill()."""
        return self.group(anchor, "O")

    # --- stereochemistry ---------------------------------------------------

    def stereocenter(self, atom: str, parity: int) -> None:
        """
        Declares configuration, ordering neighbours as currently bonded.

        Called after fill(), so hydrogens are present and the neighbour list
        is complete.
        """
        neighbors = tuple(
            neighbor for neighbor, _ in self._neighbors_of(atom)
        )
        self._centers.append(
            TetrahedralCenter(atom=atom, neighbors=neighbors, parity=parity)
        )

    def double_bond_stereo(
        self, bond_id: str, reference_source: str, reference_target: str, config: str
    ) -> None:
        self._stereo_bonds.append(
            DoubleBondStereo(
                bond=bond_id,
                reference_source=reference_source,
                reference_target=reference_target,
                config=config,
            )
        )

    # --- elaboration -------------------------------------------------------

    def _neighbors_of(self, atom_id: str) -> list[tuple[str, int]]:
        found = []
        for bond in self._bonds:
            if bond.source == atom_id:
                found.append((bond.target, bond.order))
            elif bond.target == atom_id:
                found.append((bond.source, bond.order))
        return found

    def _used(self, atom: Atom) -> int:
        bonded = sum(order for _, order in self._neighbors_of(atom.id))
        return bonded + atom.radical

    def fill(self, element: str = "H") -> None:
        """
        Saturates every unfilled valence.

        Atoms already satisfied, over-filled, or outside the oracle's model
        are left alone — over-filling is a diagnostic for the validator to
        report, not something to silently repair here.
        """
        for atom in list(self._atoms):
            if not self._oracle.supports(atom.element):
                continue

            used = self._used(atom)
            permitted = self._oracle.permitted(atom.element, atom.charge)
            target = next((v for v in permitted if v >= used), None)

            if target is None:
                continue

            for _ in range(target - used):
                self.bond(atom.id, self.atom(element))

    def freeze(self) -> Molecule:
        return Molecule(
            id=self.name.replace(" ", "-"),
            name=self.name,
            atoms=list(self._atoms),
            bonds=list(self._bonds),
            stereo_centers=list(self._centers),
            stereo_bonds=list(self._stereo_bonds),
        )


def molecule(
    name: str,
    define: Callable[[MoleculeBuilder], None],
    oracle: Optional[ValenceOracle] = None,
) -> Molecule:
    builder = MoleculeBuilder(name, oracle=oracle)
    define(builder)
    return builder.freeze()
