"""
SMILES import.

Before this, the only way to state a molecule was to write graph-construction
code, which is fine for a test corpus and hopeless for a person. `O` is water
and `c1ccccc1` is benzene, and a language nobody can write structures in is not
usable regardless of how well it checks them.

**This is a documented subset, not an authority.** RDKit is the intended
parser, exactly as it is the intended valence oracle; this exists so the
package works with no dependencies installed. Supported: the organic subset,
bracket atoms with isotope, charge, and explicit hydrogen count, bond symbols
`-` `=` `#`, branches, ring-closure digits and `%nn`, lowercase aromatics, and
`.` for disconnected components. Not supported: chirality markers `@`/`@@`
(parsed and *reported*, never silently dropped), directional bonds `/` and
`\\`, and reaction SMILES.

**Kekulisation.** The IR stores explicit alternating bond orders, so aromatic
input must be resolved to a Kekulé structure on the way in. Aromatic carbon
must take exactly one double bond; aromatic N, O, and S may take none, which
is what lets pyrrole and furan work alongside benzene. The search is
backtracking over the aromatic bonds, which is ample for rings of the size the
v0 profile targets and would need replacing for anything large.
"""

from __future__ import annotations

import re
from typing import Optional

from .elements import ValenceOracle, default_oracle
from .ir import ERROR, WARNING, Atom, Bond, Diagnostic, Molecule

ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
AROMATIC_SYMBOLS = {"b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S"}

BOND_SYMBOLS = {"-": 1, "=": 2, "#": 3, ":": 0}
AROMATIC_ORDER = 0  # placeholder, resolved by Kekulisation

# Aromatic atoms permitted to carry no double bond: their lone pair joins the
# ring system instead. Without this, pyrrole and furan are unparseable.
LONE_PAIR_AROMATIC = {"N", "O", "S"}

BRACKET_PATTERN = re.compile(
    r"^(?P<isotope>\d+)?"
    r"(?P<symbol>[A-Z][a-z]?|[bcnops])"
    r"(?P<chiral>@{1,2})?"
    r"(?:H(?P<hydrogens>\d*))?"
    r"(?P<charge>[+-]\d+|[+-]+)?$"
)


class _Parser:
    def __init__(self, text: str, oracle: ValenceOracle) -> None:
        self.text = text
        self.oracle = oracle
        self.position = 0

        self.atoms: list[Atom] = []
        self.bonds: list[tuple[int, int, int]] = []  # source, target, order
        self.aromatic: set[int] = set()
        self.explicit_hydrogens: dict[int, int] = {}
        self.diagnostics: list[Diagnostic] = []

    def fail(self, code: str, message: str) -> None:
        self.diagnostics.append(
            Diagnostic(ERROR, code, f"{message} (at position {self.position})")
        )

    def warn(self, code: str, message: str) -> None:
        self.diagnostics.append(Diagnostic(WARNING, code, message))

    # --- tokens --------------------------------------------------------

    def _add_atom(
        self,
        element: str,
        charge: int = 0,
        isotope: Optional[int] = None,
        aromatic: bool = False,
    ) -> int:
        index = len(self.atoms)
        self.atoms.append(
            Atom(
                id=f"a{index}",
                element=element,
                charge=charge,
                isotope=isotope,
            )
        )
        if aromatic:
            self.aromatic.add(index)
        return index

    def _read_bracket_atom(self) -> Optional[int]:
        closing = self.text.find("]", self.position)
        if closing == -1:
            self.fail("smiles.syntax", "unterminated '['")
            self.position = len(self.text)  # consume the rest; nothing to recover
            return None

        body = self.text[self.position + 1 : closing]
        match = BRACKET_PATTERN.match(body)

        if not match:
            self.fail("smiles.syntax", f"cannot parse bracket atom '[{body}]'")
            self.position = closing + 1
            return None

        if match.group("chiral"):
            # Reported rather than dropped: silently discarding stereochemistry
            # would produce a molecule that looks right and means something
            # else, which is the failure mode the spec cares most about.
            self.warn(
                "smiles.chiralityIgnored",
                f"chirality marker in '[{body}]' is not supported by this "
                f"parser; the resulting molecule has no configuration at that "
                f"centre",
            )

        symbol = match.group("symbol")
        aromatic = symbol in AROMATIC_SYMBOLS
        element = AROMATIC_SYMBOLS.get(symbol, symbol)

        charge = 0
        raw_charge = match.group("charge")
        if raw_charge:
            if raw_charge[-1].isdigit():
                charge = int(raw_charge[1:]) * (1 if raw_charge[0] == "+" else -1)
            else:
                charge = len(raw_charge) * (1 if raw_charge[0] == "+" else -1)

        isotope = int(match.group("isotope")) if match.group("isotope") else None

        index = self._add_atom(element, charge, isotope, aromatic)

        raw_hydrogens = match.group("hydrogens")
        if raw_hydrogens is not None:
            self.explicit_hydrogens[index] = int(raw_hydrogens or 1)
        else:
            self.explicit_hydrogens[index] = 0

        self.position = closing + 1
        return index

    def _read_organic_atom(self) -> Optional[int]:
        remaining = self.text[self.position :]

        for symbol in ("Cl", "Br"):
            if remaining.startswith(symbol):
                self.position += 2
                return self._add_atom(symbol)

        character = remaining[0]

        if character in AROMATIC_SYMBOLS:
            self.position += 1
            return self._add_atom(
                AROMATIC_SYMBOLS[character], aromatic=True
            )

        if character in ORGANIC_SUBSET:
            self.position += 1
            return self._add_atom(character)

        self.fail("smiles.syntax", f"unexpected character '{character}'")
        self.position += 1
        return None

    # --- structure -----------------------------------------------------

    def parse(self) -> None:
        previous: Optional[int] = None
        pending_bond: Optional[int] = None
        branches: list[Optional[int]] = []
        ring_bonds: dict[str, tuple[int, Optional[int]]] = {}

        last_position = -1

        while self.position < len(self.text):
            # Every branch below must consume at least one character. A branch
            # that reports an error and forgets to advance turns a malformed
            # input into a hang, so the invariant is enforced rather than
            # trusted — an unterminated '[' did exactly that.
            if self.position == last_position:
                self.fail(
                    "smiles.internal",
                    "parser made no progress; refusing to loop",
                )
                break
            last_position = self.position

            character = self.text[self.position]

            if character in "/\\":
                self.warn(
                    "smiles.directionalBondIgnored",
                    "directional bonds are not supported; E/Z configuration "
                    "is not recovered from this input",
                )
                self.position += 1
                continue

            if character == "(":
                branches.append(previous)
                self.position += 1
                continue

            if character == ")":
                if not branches:
                    self.fail("smiles.syntax", "unmatched ')'")
                    self.position += 1
                    continue
                previous = branches.pop()
                self.position += 1
                continue

            if character in BOND_SYMBOLS:
                pending_bond = BOND_SYMBOLS[character]
                self.position += 1
                continue

            if character == ".":
                previous = None
                pending_bond = None
                self.position += 1
                continue

            if character.isdigit() or character == "%":
                label = self._read_ring_label()
                if label is None:
                    continue

                if label in ring_bonds:
                    partner, stored_bond = ring_bonds.pop(label)
                    order = pending_bond or stored_bond
                    if order is None:
                        order = (
                            AROMATIC_ORDER
                            if partner in self.aromatic and previous in self.aromatic
                            else 1
                        )
                    if previous is not None:
                        self.bonds.append((partner, previous, order))
                else:
                    if previous is None:
                        self.fail("smiles.syntax", "ring closure before any atom")
                    else:
                        ring_bonds[label] = (previous, pending_bond)

                pending_bond = None
                continue

            index = (
                self._read_bracket_atom()
                if character == "["
                else self._read_organic_atom()
            )

            if index is None:
                continue

            if previous is not None:
                order = pending_bond
                if order is None:
                    order = (
                        AROMATIC_ORDER
                        if previous in self.aromatic and index in self.aromatic
                        else 1
                    )
                self.bonds.append((previous, index, order))

            previous = index
            pending_bond = None

        if ring_bonds:
            self.fail(
                "smiles.syntax",
                f"unclosed ring bond(s): {sorted(ring_bonds)}",
            )
        if branches:
            self.fail("smiles.syntax", "unmatched '('")

    def _read_ring_label(self) -> Optional[str]:
        if self.text[self.position] == "%":
            label = self.text[self.position + 1 : self.position + 3]
            if len(label) != 2 or not label.isdigit():
                self.fail("smiles.syntax", "'%' must be followed by two digits")
                self.position += 1
                return None
            self.position += 3
            return label

        label = self.text[self.position]
        self.position += 1
        return label

    # --- kekulisation --------------------------------------------------

    def kekulise(self) -> bool:
        """
        Assigns alternating orders across the aromatic system.

        Aromatic carbon must end with exactly one double bond; N, O, and S may
        end with none. Returns False when no assignment exists, which is a
        refusal rather than a guess.
        """
        aromatic_indices = [
            position
            for position, (_, _, order) in enumerate(self.bonds)
            if order == AROMATIC_ORDER
        ]

        if not aromatic_indices:
            return True

        incident: dict[int, list[int]] = {}
        for position in aromatic_indices:
            source, target, _ = self.bonds[position]
            incident.setdefault(source, []).append(position)
            incident.setdefault(target, []).append(position)

        doubled: dict[int, bool] = {}

        def satisfied(atom_index: int) -> bool:
            count = sum(
                1 for position in incident[atom_index] if doubled.get(position)
            )
            if count > 1:
                return False
            if count == 1:
                return True
            element = self.atoms[atom_index].element
            return element in LONE_PAIR_AROMATIC

        def search(cursor: int) -> bool:
            if cursor == len(aromatic_indices):
                return all(satisfied(atom) for atom in incident)

            position = aromatic_indices[cursor]
            source, target, _ = self.bonds[position]

            for choice in (True, False):
                doubled[position] = choice

                # Prune: an atom whose bonds are all decided must be valid now.
                decided = all(
                    other in doubled
                    for atom in (source, target)
                    for other in incident[atom]
                )
                if not decided or (satisfied(source) and satisfied(target)):
                    if search(cursor + 1):
                        return True

            del doubled[position]
            return False

        if not search(0):
            return False

        for position in aromatic_indices:
            source, target, _ = self.bonds[position]
            self.bonds[position] = (
                source,
                target,
                2 if doubled.get(position) else 1,
            )

        return True

    # --- elaboration ---------------------------------------------------

    def build(self, name: str) -> Molecule:
        """Adds implicit hydrogens and freezes the graph."""
        atoms = list(self.atoms)
        bonds = [
            Bond(id=f"b{position}", source=f"a{s}", target=f"a{t}", order=order)
            for position, (s, t, order) in enumerate(self.bonds)
        ]

        used: dict[int, int] = {index: 0 for index in range(len(atoms))}
        for source, target, order in self.bonds:
            used[source] += order
            used[target] += order

        counter = len(bonds)
        for index, atom in enumerate(list(atoms)):
            if index in self.explicit_hydrogens:
                wanted = self.explicit_hydrogens[index]
            elif self.oracle.supports(atom.element):
                permitted = self.oracle.permitted(atom.element, atom.charge)
                target_valence = next(
                    (value for value in permitted if value >= used[index]), None
                )
                wanted = (
                    0 if target_valence is None else target_valence - used[index]
                )
            else:
                wanted = 0

            for _ in range(wanted):
                hydrogen_index = len(atoms)
                atoms.append(Atom(id=f"a{hydrogen_index}", element="H"))
                bonds.append(
                    Bond(
                        id=f"b{counter}",
                        source=atom.id,
                        target=f"a{hydrogen_index}",
                        order=1,
                    )
                )
                counter += 1

        return Molecule(id=name.replace(" ", "-"), name=name, atoms=atoms, bonds=bonds)


def parse(
    text: str,
    name: Optional[str] = None,
    oracle: Optional[ValenceOracle] = None,
) -> tuple[Optional[Molecule], list[Diagnostic]]:
    """
    Reads a SMILES string.

    Returns the molecule alongside anything noteworthy about the reading —
    including warnings for stereochemistry this parser cannot represent, so a
    caller is never silently handed a structure that means less than its input
    did.
    """
    if not text or not text.strip():
        return None, [
            Diagnostic(ERROR, "smiles.syntax", "empty SMILES string")
        ]

    parser = _Parser(text.strip(), oracle or default_oracle())
    parser.parse()

    if any(d.severity == ERROR for d in parser.diagnostics):
        return None, parser.diagnostics

    if not parser.kekulise():
        parser.diagnostics.append(
            Diagnostic(
                ERROR,
                "chemistry.unsupported",
                "cannot resolve the aromatic system to alternating bond "
                "orders; supply an explicit Kekulé structure instead",
            )
        )
        return None, parser.diagnostics

    return parser.build(name or text.strip()), parser.diagnostics


def parse_or_raise(
    text: str,
    name: Optional[str] = None,
    oracle: Optional[ValenceOracle] = None,
) -> Molecule:
    """Convenience for known-good input, as in tests and demos."""
    molecule, diagnostics = parse(text, name, oracle)
    if molecule is None:
        raise ValueError(
            f"[pouring] cannot parse SMILES {text!r}: "
            f"{'; '.join(d.message for d in diagnostics)}"
        )
    return molecule
