"""
Interchange: reading the IR back, and reading structures in.

Two claims are under test here. That the canonical JSON encoding is genuinely
an interchange format rather than a one-way dump — spec 13.1 lists round trip
as a conformance class, and until a decoder existed it could not be run. And
that a structure written the way a chemist writes one arrives at the same
molecule as the hand-built graph, which is the cross-check that two
independent construction paths agree.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pouring import (  # noqa: E402
    CRN,
    SYNTHESIS,
    Document,
    Route,
    approximate_majority,
    aspirin,
    aspirin_synthesis,
    autocatalysis,
    benzene,
    butene,
    carbon_dioxide,
    check_structure,
    combustion,
    content_hash,
    default_oracle,
    dumps,
    formula,
    glucose,
    hydrogen_peroxide,
    loads,
    methane,
    molar_mass,
    parse_smiles,
    ring_count,
    round_trip,
    water,
)
from pouring.smiles import parse_or_raise  # noqa: E402

ORACLE = default_oracle()


def wrap(body):
    profile = CRN if body.kind == "reactionNetwork" else SYNTHESIS
    return Document(target_profile=profile, body=body)


def as_route(mol, name="probe"):
    return Route(
        id=name,
        name=name,
        target=mol.id,
        starting_materials=[mol.id],
        species=[mol],
        reactions=[],
    )


class TestRoundTrip(unittest.TestCase):
    def test_bodies_survive_encoding(self):
        for body in [
            autocatalysis(),
            approximate_majority(10, 5),
            aspirin_synthesis(),
            combustion(),
        ]:
            document = wrap(body)
            self.assertEqual(
                dumps(document), dumps(round_trip(document)), body.name
            )

    def test_structure_survives(self):
        document = wrap(as_route(glucose(), "glucose-probe"))
        recovered = round_trip(document).body.species[0]
        self.assertEqual(formula(recovered), "C6H12O6")
        self.assertEqual(molar_mass(recovered), molar_mass(glucose()))

    def test_stereochemistry_survives(self):
        """
        The point of encoding stereo at all: an enantiomer that loses its
        configuration in transit becomes a different claim than the one sent.
        """
        document = wrap(as_route(glucose(), "glucose-probe"))
        recovered = round_trip(document).body.species[0]
        self.assertEqual(len(recovered.stereo_centers), 4)
        self.assertEqual(content_hash(recovered), content_hash(glucose()))
        self.assertNotEqual(
            content_hash(recovered), content_hash(glucose(mirrored=True))
        )

    def test_geometric_configuration_survives(self):
        document = wrap(as_route(butene(), "butene-probe"))
        recovered = round_trip(document).body.species[0]
        self.assertEqual(content_hash(recovered), content_hash(butene()))

    def test_conditions_and_yield_survive(self):
        recovered = round_trip(wrap(aspirin_synthesis())).body.reactions[0]
        self.assertEqual(recovered.conditions.temperature.unit, "degC")
        self.assertEqual(recovered.conditions.catalyst, "sulfuric acid")
        self.assertEqual(recovered.reaction_yield.product, "aspirin")

    def test_exact_coefficients_survive(self):
        recovered = round_trip(wrap(approximate_majority(10, 5))).body
        self.assertEqual(recovered.reactions[0].products[0].coefficient, 2)


class TestDecoderRejectsBadInput(unittest.TestCase):
    """Malformed input is ordinary, so it produces diagnostics, not tracebacks."""

    def test_invalid_json(self):
        document, diagnostics = loads("not json at all")
        self.assertIsNone(document)
        self.assertEqual(diagnostics[0].code, "schema.malformed")

    def test_unsupported_schema_version(self):
        document, diagnostics = loads('{"schemaVersion": "9.9.9"}')
        self.assertIsNone(document)
        self.assertEqual(diagnostics[0].code, "schema.version")

    def test_unknown_target_profile(self):
        document, diagnostics = loads(
            '{"schemaVersion": "0.1.0", "targetProfile": "nonsense"}'
        )
        self.assertIsNone(document)
        self.assertEqual(diagnostics[0].code, "schema.version")

    def test_unknown_body_kind(self):
        document, diagnostics = loads(
            '{"schemaVersion": "0.1.0", "targetProfile": "pouring:crn-v0",'
            ' "body": {"kind": "sculpture"}}'
        )
        self.assertIsNone(document)
        self.assertEqual(diagnostics[0].code, "schema.unknownKind")

    def test_missing_required_field(self):
        document, diagnostics = loads(
            '{"schemaVersion": "0.1.0", "targetProfile": "pouring:crn-v0",'
            ' "body": {"kind": "reactionNetwork", "id": "x",'
            ' "species": [{"kind": "molecule", "name": "no id"}]}}'
        )
        self.assertIn("schema.malformed", {d.code for d in diagnostics})


class TestSmiles(unittest.TestCase):
    def test_organic_subset(self):
        for text, expected in [
            ("O", "H2O"),
            ("OO", "H2O2"),
            ("C", "CH4"),
            ("CCO", "C2H6O"),
            ("O=C=O", "CO2"),
            ("C#N", "CHN"),
            ("CC(C)C", "C4H10"),
        ]:
            self.assertEqual(formula(parse_or_raise(text)), expected, text)

    def test_rings(self):
        self.assertEqual(formula(parse_or_raise("C1CCCCC1")), "C6H12")
        self.assertEqual(ring_count(parse_or_raise("C1CCCCC1")), 1)

    def test_aromatics_are_kekulised(self):
        """
        The IR stores explicit alternating orders, so lowercase aromatic input
        must be resolved on the way in or refused.
        """
        for text, expected in [
            ("c1ccccc1", "C6H6"),
            ("c1ccncc1", "C5H5N"),
            ("c1cc[nH]c1", "C4H5N"),
            ("c1ccoc1", "C4H4O"),
        ]:
            molecule = parse_or_raise(text)
            self.assertEqual(formula(molecule), expected, text)
            self.assertEqual(
                [d for d in check_structure(molecule, ORACLE) if d.severity == "error"],
                [],
                text,
            )

    def test_bracket_atoms(self):
        self.assertEqual(parse_or_raise("[Na+]").atoms[0].charge, 1)
        self.assertEqual(parse_or_raise("[O-]C").atoms[0].charge, -1)
        self.assertEqual(formula(parse_or_raise("[NH4+]")), "H4N")
        self.assertEqual(parse_or_raise("[13CH4]").atoms[0].isotope, 13)

    def test_isotope_changes_mass(self):
        self.assertGreater(
            molar_mass(parse_or_raise("[13CH4]")), molar_mass(parse_or_raise("C"))
        )

    def test_agrees_with_hand_built_molecules(self):
        """
        Two independent construction paths must converge on one identity. If
        they ever disagree, one of the parser and the builder is wrong and the
        canonical form cannot arbitrate.
        """
        for text, built in [
            ("O", water()),
            ("OO", hydrogen_peroxide()),
            ("C", methane()),
            ("O=C=O", carbon_dioxide()),
            ("c1ccccc1", benzene()),
            ("C1=CC=CC=C1", benzene()),
            ("CC(=O)Oc1ccccc1C(=O)O", aspirin()),
        ]:
            self.assertEqual(
                content_hash(parse_or_raise(text)), content_hash(built), text
            )

    def test_both_kekule_inputs_agree(self):
        self.assertEqual(
            content_hash(parse_or_raise("C1=CC=CC=C1")),
            content_hash(parse_or_raise("c1ccccc1")),
        )

    def test_parsed_molecules_survive_round_trip(self):
        molecule = parse_or_raise("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
        recovered = round_trip(wrap(as_route(molecule))).body.species[0]
        self.assertEqual(content_hash(recovered), content_hash(molecule))

    def test_syntax_errors_are_reported(self):
        for text in ["C1CC", "CC)", "[C", "", "CC(C"]:
            molecule, diagnostics = parse_smiles(text)
            self.assertIsNone(molecule, text)
            self.assertTrue(diagnostics, text)

    def test_chirality_is_reported_not_dropped(self):
        """
        Silently discarding a stereo marker would hand back a molecule that
        looks right and means something else — precisely the failure the spec
        treats as disqualifying.
        """
        molecule, diagnostics = parse_smiles("N[C@@H](C)C(=O)O")
        self.assertIsNotNone(molecule)
        self.assertIn(
            "smiles.chiralityIgnored", {d.code for d in diagnostics}
        )
        self.assertEqual(molecule.stereo_centers, [])

    def test_disconnected_components_parse_but_do_not_validate(self):
        molecule, _ = parse_smiles("[Na+].[Cl-]")
        self.assertIsNotNone(molecule)
        self.assertIn(
            "molecule.disconnected",
            {d.code for d in check_structure(molecule, ORACLE)},
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
