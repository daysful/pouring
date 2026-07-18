"""
Reading and writing the canonical JSON encoding (spec 10.1).

The IR claims to be language-neutral data that many frontends produce and many
backends consume. Until this module existed that claim was half true: documents
could be written and nothing could read them back, which makes an interchange
format that does not interchange.

Decoding reports problems as diagnostics rather than raising. A malformed
document is an ordinary input to a compiler, not an exceptional condition, and
a caller loading someone else's file deserves a list of what is wrong with it
rather than a stack trace at the first bad field.
"""

from __future__ import annotations

import json
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from typing import Any, Optional

from .ir import (
    BALANCED,
    CHEMISTRY_V0,
    CRN,
    ERROR,
    NORMALIZE_V0,
    SCHEMA_VERSION,
    SYNTHESIS,
    Atom,
    Bond,
    Component,
    Conditions,
    Diagnostic,
    Document,
    DoubleBondStereo,
    FormalSpecies,
    Kinetics,
    Molecule,
    MultiComponentSpecies,
    Participant,
    Reaction,
    ReactionNetwork,
    Route,
    TetrahedralCenter,
    Yield,
)
from .units import Quantity, UnknownUnit

SUPPORTED_SCHEMA_VERSIONS = {SCHEMA_VERSION}


class _Reader:
    """Accumulates diagnostics while walking a document."""

    def __init__(self) -> None:
        self.diagnostics: list[Diagnostic] = []

    def fail(self, code: str, message: str, path: str | None = None) -> None:
        self.diagnostics.append(Diagnostic(ERROR, code, message, path=path))

    def require(self, data: Any, key: str, path: str) -> Any:
        if not isinstance(data, dict) or key not in data:
            self.fail(
                "schema.malformed", f"missing required field '{key}'", path
            )
            return None
        return data[key]

    @property
    def ok(self) -> bool:
        return not self.diagnostics


# --- scalars ---------------------------------------------------------------


def _rational(data: Any, reader: _Reader, path: str) -> Fraction:
    """
    Spec 5.1: rationals travel as string numerator and denominator, so an
    exact coefficient survives the trip. A bare number is accepted as a
    convenience for hand-written documents.
    """
    if isinstance(data, (int, float)):
        return Fraction(str(data))

    if isinstance(data, dict):
        try:
            return Fraction(
                int(data.get("numerator", 1)), int(data.get("denominator", 1))
            )
        except (ValueError, ZeroDivisionError):
            reader.fail("schema.malformed", f"invalid rational {data!r}", path)
            return Fraction(1)

    reader.fail("schema.malformed", f"expected a rational, got {data!r}", path)
    return Fraction(1)


def _quantity(data: Any, reader: _Reader, path: str) -> Optional[Quantity]:
    if not isinstance(data, dict):
        reader.fail("schema.malformed", f"expected a quantity, got {data!r}", path)
        return None

    unit = data.get("unit")
    if not isinstance(unit, str):
        reader.fail("schema.malformed", "quantity has no unit", path)
        return None

    try:
        value = Decimal(str(data.get("value", "0")))
    except InvalidOperation:
        reader.fail(
            "schema.malformed", f"invalid decimal {data.get('value')!r}", path
        )
        return None

    try:
        return Quantity(
            value=value,
            unit=unit,
            uncertainty=(
                Decimal(str(data["uncertainty"]))
                if data.get("uncertainty") is not None
                else None
            ),
            qualifier=data.get("qualifier", "exact"),
        )
    except UnknownUnit as error:
        reader.fail("unit.unknown", str(error), path)
        return None


# --- structure -------------------------------------------------------------


def _atom(data: dict, reader: _Reader, path: str) -> Optional[Atom]:
    atom_id = reader.require(data, "id", path)
    element = reader.require(data, "element", path)
    if atom_id is None or element is None:
        return None

    return Atom(
        id=str(atom_id),
        element=str(element),
        charge=int(data.get("charge", 0)),
        radical=int(data.get("radical", 0)),
        isotope=data.get("isotope"),
    )


def _bond(data: dict, reader: _Reader, path: str) -> Optional[Bond]:
    bond_id = reader.require(data, "id", path)
    source = reader.require(data, "from", path)
    target = reader.require(data, "to", path)
    if bond_id is None or source is None or target is None:
        return None

    return Bond(
        id=str(bond_id),
        source=str(source),
        target=str(target),
        order=int(data.get("order", 1)),
    )


def _molecule(data: dict, reader: _Reader) -> Optional[Molecule]:
    identifier = reader.require(data, "id", "molecule")
    if identifier is None:
        return None

    path = str(identifier)

    atoms = [
        atom
        for atom in (
            _atom(entry, reader, path) for entry in data.get("atoms", [])
        )
        if atom is not None
    ]
    bonds = [
        bond
        for bond in (
            _bond(entry, reader, path) for entry in data.get("bonds", [])
        )
        if bond is not None
    ]

    centers = [
        TetrahedralCenter(
            atom=str(entry["atom"]),
            neighbors=tuple(str(n) for n in entry.get("neighbors", ())),
            parity=int(entry.get("parity", 1)),
        )
        for entry in data.get("stereoCenters", [])
        if "atom" in entry
    ]

    stereo_bonds = [
        DoubleBondStereo(
            bond=str(entry["bond"]),
            reference_source=str(entry.get("referenceFrom", "")),
            reference_target=str(entry.get("referenceTo", "")),
            config=str(entry.get("config", "cis")),
        )
        for entry in data.get("stereoBonds", [])
        if "bond" in entry
    ]

    return Molecule(
        id=path,
        name=str(data.get("name", path)),
        atoms=atoms,
        bonds=bonds,
        stereo_centers=centers,
        stereo_bonds=stereo_bonds,
    )


def _species(data: Any, reader: _Reader):
    if not isinstance(data, dict):
        reader.fail("schema.malformed", f"expected a species, got {data!r}")
        return None

    kind = data.get("kind")

    if kind == "molecule":
        return _molecule(data, reader)

    if kind == "formalSpecies":
        identifier = reader.require(data, "id", "formalSpecies")
        if identifier is None:
            return None
        return FormalSpecies(
            id=str(identifier),
            name=str(data.get("name", identifier)),
            binding=data.get("binding"),
        )

    if kind == "multiComponentSpecies":
        identifier = reader.require(data, "id", "multiComponentSpecies")
        if identifier is None:
            return None
        return MultiComponentSpecies(
            id=str(identifier),
            name=str(data.get("name", identifier)),
            components=[
                Component(
                    entity=str(entry.get("entity", "")),
                    ratio=_rational(entry.get("ratio"), reader, str(identifier)),
                )
                for entry in data.get("components", [])
            ],
        )

    reader.fail("schema.unknownKind", f"unknown species kind {kind!r}")
    return None


# --- transformation --------------------------------------------------------


def _participant(data: dict, reader: _Reader, path: str) -> Participant:
    return Participant(
        species=str(data.get("species", "")),
        coefficient=_rational(data.get("coefficient"), reader, path),
    )


def _conditions(data: Any, reader: _Reader, path: str) -> Optional[Conditions]:
    if not isinstance(data, dict):
        return None

    return Conditions(
        temperature=(
            _quantity(data["temperature"], reader, path)
            if "temperature" in data
            else None
        ),
        time=_quantity(data["time"], reader, path) if "time" in data else None,
        pressure=(
            _quantity(data["pressure"], reader, path) if "pressure" in data else None
        ),
        solvent=data.get("solvent"),
        catalyst=data.get("catalyst"),
    )


def _reaction(data: dict, reader: _Reader) -> Optional[Reaction]:
    identifier = reader.require(data, "id", "reaction")
    if identifier is None:
        return None

    path = str(identifier)

    kinetics = None
    if isinstance(data.get("kinetics"), dict):
        kinetics = Kinetics(
            rate_constant=float(data["kinetics"].get("rateConstant", 1.0)),
            rate_law=str(data["kinetics"].get("rateLaw", "massAction")),
        )

    recorded_yield = None
    if isinstance(data.get("yield"), dict):
        entry = data["yield"]
        try:
            recorded_yield = Yield(
                product=str(entry.get("product", "")),
                basis=str(entry.get("basis", "")),
                value=Decimal(str(entry.get("value", "0"))),
                kind=str(entry.get("kind", "isolated")),
                method=entry.get("method"),
            )
        except InvalidOperation:
            reader.fail("schema.malformed", "invalid yield value", path)

    return Reaction(
        id=path,
        name=str(data.get("name", path)),
        reactants=[
            _participant(entry, reader, path) for entry in data.get("reactants", [])
        ],
        products=[
            _participant(entry, reader, path) for entry in data.get("products", [])
        ],
        kinetics=kinetics,
        conditions=_conditions(data.get("conditions"), reader, path),
        equation_status=str(data.get("equationStatus", BALANCED)),
        reaction_yield=recorded_yield,
        atom_map=data.get("atomMap"),
    )


def _body(data: Any, reader: _Reader):
    if not isinstance(data, dict):
        reader.fail("schema.malformed", "document has no body")
        return None

    kind = data.get("kind")

    species = [
        entry
        for entry in (_species(item, reader) for item in data.get("species", []))
        if entry is not None
    ]
    reactions = [
        entry
        for entry in (_reaction(item, reader) for item in data.get("reactions", []))
        if entry is not None
    ]

    identifier = str(data.get("id", "body"))

    if kind == "reactionNetwork":
        return ReactionNetwork(
            id=identifier,
            name=str(data.get("name", identifier)),
            species=species,
            reactions=reactions,
            inputs={
                str(key): float(value)
                for key, value in (data.get("inputs") or {}).items()
            },
            outputs=[str(entry) for entry in data.get("outputs", [])],
            semantics=data.get("semantics"),
        )

    if kind == "route":
        return Route(
            id=identifier,
            name=str(data.get("name", identifier)),
            target=str(data.get("target", "")),
            starting_materials=[
                str(entry) for entry in data.get("startingMaterials", [])
            ],
            species=species,
            reactions=reactions,
        )

    reader.fail("schema.unknownKind", f"unknown body kind {kind!r}")
    return None


# --- documents -------------------------------------------------------------


def decode_document(data: Any) -> tuple[Optional[Document], list[Diagnostic]]:
    """
    Reads a document, returning it alongside whatever was wrong with it.

    A `None` document means the input could not be understood at all; a
    document with diagnostics means it was understood and is defective, which
    are different situations and should not be conflated.
    """
    reader = _Reader()

    if not isinstance(data, dict):
        reader.fail("schema.malformed", "document is not an object")
        return None, reader.diagnostics

    version = str(data.get("schemaVersion", ""))
    if version not in SUPPORTED_SCHEMA_VERSIONS:
        reader.fail(
            "schema.version",
            f"unsupported schemaVersion {version!r}; this build reads "
            f"{sorted(SUPPORTED_SCHEMA_VERSIONS)}",
        )
        return None, reader.diagnostics

    target = str(data.get("targetProfile", ""))
    if target not in (SYNTHESIS, CRN):
        reader.fail("schema.version", f"unknown targetProfile {target!r}")
        return None, reader.diagnostics

    body = _body(data.get("body"), reader)
    if body is None:
        return None, reader.diagnostics

    document = Document(
        target_profile=target,
        body=body,
        chemistry_profile=str(data.get("chemistryProfile", CHEMISTRY_V0)),
        normalization_profile=str(
            data.get("normalizationProfile", NORMALIZE_V0)
        ),
        schema_version=version,
    )

    return document, reader.diagnostics


def encode_document(document: Document) -> dict:
    return document.to_json()


def loads(text: str) -> tuple[Optional[Document], list[Diagnostic]]:
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        return None, [
            Diagnostic(ERROR, "schema.malformed", f"invalid JSON: {error}")
        ]
    return decode_document(data)


def dumps(document: Document, indent: int | None = 2) -> str:
    return json.dumps(encode_document(document), indent=indent, sort_keys=True)


def round_trip(document: Document) -> Document:
    """
    Encodes and decodes a document.

    Useful as a property: `dumps(round_trip(d)) == dumps(d)` must hold for any
    document the implementation can produce, or the encoding is lossy.
    """
    decoded, diagnostics = loads(dumps(document))
    if decoded is None:
        raise ValueError(
            f"[pouring] round trip failed: "
            f"{'; '.join(d.message for d in diagnostics)}"
        )
    return decoded
