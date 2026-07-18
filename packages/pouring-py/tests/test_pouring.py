"""
Conformance tests.

Organised around the classes in spec 13.1 — golden cases, permutation and
bond-reversal invariance, round trips, canonicalisation, stereo distinction,
and unit metamorphics — plus the cross-target case that keeps "one shared IR"
honest.

    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pouring import (  # noqa: E402
    CRN,
    SYNTHESIS,
    Atom,
    Bond,
    Component,
    Document,
    MultiComponentSpecies,
    ProfileValenceOracle,
    Reaction,
    Route,
    acetic_acid,
    approximate_majority,
    approximate_majority_as_route,
    aspirin,
    aspirin_synthesis,
    autocatalysis,
    balance,
    balanced_equation,
    benzene,
    butene,
    canonical_encoding,
    carbon_dioxide,
    check_stereo,
    check_structure,
    combustion,
    connected_components,
    content_hash,
    convert,
    default_oracle,
    describe,
    dioxygen,
    errors,
    formula,
    glucose,
    heavy_water,
    hydrogen_peroxide,
    hydroperoxyl,
    is_valid,
    majority_verdict,
    methane,
    mirror,
    molar_mass,
    molecule,
    perceive_aromatic_bonds,
    potential_stereocenters,
    quantity,
    ring_count,
    salicylic_acid,
    simulate_deterministic,
    simulate_stochastic,
    sodium_chloride,
    validate,
    warnings,
    water,
)
from pouring.units import DimensionMismatch, UnknownUnit, add  # noqa: E402

ORACLE = ProfileValenceOracle()


def codes(diagnostics) -> set[str]:
    return {d.code for d in diagnostics}


def structural_errors(mol):
    return [d for d in check_structure(mol, ORACLE) if d.severity == "error"]


# --- the shared-IR claim ---------------------------------------------------


class TestCrossTarget(unittest.TestCase):
    """One reaction set, two readings, exactly three differences."""

    def setUp(self):
        self.as_network = validate(
            Document(target_profile=CRN, body=approximate_majority(10, 5))
        )
        self.as_route = validate(
            Document(
                target_profile=SYNTHESIS, body=approximate_majority_as_route()
            )
        )

    def test_valid_as_a_program(self):
        self.assertTrue(is_valid(self.as_network))

    def test_invalid_as_a_plan(self):
        self.assertFalse(is_valid(self.as_route))

    def test_formal_species_normal_in_one_unresolved_in_other(self):
        self.assertIn("species.formal", codes(self.as_network))
        self.assertIn("species.unresolved", codes(self.as_route))
        self.assertNotIn("species.unresolved", codes(self.as_network))

    def test_cycles_are_the_mechanism_not_an_error(self):
        self.assertNotIn("route.cycle", codes(self.as_network))
        self.assertIn("route.cycle", codes(self.as_route))

    def test_no_fourth_difference(self):
        """
        Codes appearing on one side and not the other must all be licensed by
        spec 3.1. A new divergence means the readings have drifted apart in a
        way the shared-IR claim does not allow.
        """
        licensed = {
            "species.formal",
            "species.unresolved",
            "route.cycle",
            "route.unusedSpecies",
            "route.unconsumedIntermediate",
            "kinetics.missing",
            "kinetics.semanticsUnspecified",
            "reaction.unbalanced",
            "reaction.underdetermined",
            "reaction.charge",
        }
        difference = codes(self.as_network) ^ codes(self.as_route)
        self.assertTrue(
            difference <= licensed,
            f"unlicensed divergence: {difference - licensed}",
        )


class TestBalanceIsProfileBound(unittest.TestCase):
    def test_autocatalysis_permitted_under_crn(self):
        result = validate(Document(target_profile=CRN, body=autocatalysis()))
        self.assertNotIn("reaction.unbalanced", codes(result))
        self.assertTrue(is_valid(result))

    def test_real_synthesis_balances(self):
        result = validate(
            Document(target_profile=SYNTHESIS, body=aspirin_synthesis())
        )
        self.assertTrue(is_valid(result), [str(d) for d in errors(result)])

    def test_combustion_balances(self):
        result = validate(Document(target_profile=SYNTHESIS, body=combustion()))
        self.assertTrue(is_valid(result), [str(d) for d in errors(result)])

    def test_checker_actually_fires(self):
        result = validate(
            Document(
                target_profile=SYNTHESIS, body=aspirin_synthesis(sabotage=True)
            )
        )
        self.assertIn("reaction.unbalanced", codes(result))


class TestKineticsAreProfileBound(unittest.TestCase):
    def _stripped_network(self):
        network = approximate_majority(10, 5)
        for reaction in network.reactions:
            reaction.kinetics = None
        network.semantics = None
        return network

    def _stripped_route(self):
        route = approximate_majority_as_route()
        for reaction in route.reactions:
            reaction.kinetics = None
        return route

    def test_required_under_crn(self):
        result = validate(
            Document(target_profile=CRN, body=self._stripped_network())
        )
        self.assertIn("kinetics.missing", codes(result))
        self.assertIn("kinetics.semanticsUnspecified", codes(result))

    def test_optional_under_synthesis(self):
        result = validate(
            Document(target_profile=SYNTHESIS, body=self._stripped_route())
        )
        self.assertNotIn("kinetics.missing", codes(result))


# --- structure -------------------------------------------------------------


class TestStructure(unittest.TestCase):
    def test_formulas(self):
        expected = {
            "H2O": water(),
            "H2O2": hydrogen_peroxide(),
            "HO2": hydroperoxyl(),
            "C6H6": benzene(),
            "C6H12O6": glucose(),
            "C7H6O3": salicylic_acid(),
            "C9H8O4": aspirin(),
            "C2H4O2": acetic_acid(),
            "C4H8": butene(),
            "CH4": methane(),
            "CO2": carbon_dioxide(),
        }
        for wanted, mol in expected.items():
            self.assertEqual(formula(mol), wanted, mol.name)

    def test_masses(self):
        self.assertAlmostEqual(float(molar_mass(water())), 18.015, places=3)
        self.assertAlmostEqual(float(molar_mass(benzene())), 78.114, places=3)
        self.assertAlmostEqual(float(molar_mass(glucose())), 180.156, places=3)

    def test_isotopes_change_mass(self):
        """D2O is 20.03 g/mol, not water's 18.02."""
        self.assertAlmostEqual(float(molar_mass(heavy_water())), 20.028, places=2)
        self.assertGreater(molar_mass(heavy_water()), molar_mass(water()))

    def test_benzene_has_a_ring(self):
        self.assertEqual(ring_count(benzene()), 1)
        self.assertEqual(ring_count(water()), 0)
        self.assertEqual(ring_count(aspirin()), 1)

    def test_corpus_is_structurally_valid(self):
        for mol in [
            water(),
            hydrogen_peroxide(),
            hydroperoxyl(),
            benzene(),
            benzene(alternate=True),
            glucose(),
            glucose(mirrored=True),
            salicylic_acid(),
            aspirin(),
            acetic_acid(),
            butene(),
            methane(),
            dioxygen(),
            carbon_dioxide(),
        ]:
            self.assertEqual(structural_errors(mol), [], describe(mol))

    def test_connectivity(self):
        self.assertEqual(len(connected_components(water())), 1)

    def test_disconnected_is_rejected(self):
        def define(m):
            m.atom("O")
            m.atom("O")  # two atoms, no bond
            m.fill()

        found = check_structure(molecule("split", define), ORACLE)
        self.assertIn("molecule.disconnected", codes(found))


class TestValence(unittest.TestCase):
    """Spec 4.5 — the rule, and its explicit limits."""

    def test_floating_valence_is_an_error(self):
        def define(m):
            oxygen = m.atom("O")
            m.bond(oxygen, m.atom("H"))  # only one H: O is short by one

        found = check_structure(molecule("OH fragment", define), ORACLE)
        self.assertIn("valence.floating", codes(found))

    def test_radical_must_be_declared(self):
        """
        Hydroperoxyl validates only because its unpaired electron is stated.
        Remove the declaration and the same graph is an error, not a silent
        radical.
        """
        self.assertEqual(structural_errors(hydroperoxyl()), [])

        undeclared = hydroperoxyl()
        undeclared.atoms = [
            Atom(id=a.id, element=a.element, charge=a.charge, radical=0)
            for a in undeclared.atoms
        ]
        self.assertIn("valence.floating", codes(check_structure(undeclared, ORACLE)))

    def test_exceeded_valence_is_an_error(self):
        def define(m):
            carbon = m.atom("C")
            for _ in range(5):
                m.bond(carbon, m.atom("H"))

        found = check_structure(molecule("CH5", define), ORACLE)
        self.assertIn("valence.exceeded", codes(found))

    def test_charge_shifts_valence_by_element_type(self):
        """
        Ammonium nitrogen takes four bonds; both carbocation and carbanion
        take three. The direction of the shift depends on the element.
        """
        self.assertIn(4, ORACLE.permitted("N", 1))
        self.assertIn(1, ORACLE.permitted("O", -1))
        self.assertIn(3, ORACLE.permitted("C", 1))
        self.assertIn(3, ORACLE.permitted("C", -1))

    def test_self_loop_rejected(self):
        mol = water()
        mol.bonds = list(mol.bonds) + [
            Bond(id="bad", source=mol.atoms[0].id, target=mol.atoms[0].id)
        ]
        self.assertIn("bond.selfLoop", codes(check_structure(mol, ORACLE)))

    def test_duplicate_bond_rejected(self):
        mol = water()
        first = mol.bonds[0]
        mol.bonds = list(mol.bonds) + [
            Bond(id="dup", source=first.source, target=first.target)
        ]
        self.assertIn("bond.duplicate", codes(check_structure(mol, ORACLE)))


# --- canonicalisation ------------------------------------------------------


class TestCanonical(unittest.TestCase):
    def test_kekule_forms_are_the_same_molecule(self):
        """Spec 4.2: both drawings of benzene must canonicalise identically."""
        self.assertEqual(
            content_hash(benzene()), content_hash(benzene(alternate=True))
        )

    def test_aromatic_ring_is_perceived(self):
        self.assertEqual(len(perceive_aromatic_bonds(benzene())), 6)
        self.assertEqual(len(perceive_aromatic_bonds(water())), 0)

    def test_permutation_invariance(self):
        """Atom declaration order must not change identity."""

        def first(m):
            oxygen = m.atom("O")
            m.bond(oxygen, m.atom("H"))
            m.bond(oxygen, m.atom("H"))

        def second(m):
            hydrogen = m.atom("H")
            oxygen = m.atom("O")
            other = m.atom("H")
            m.bond(oxygen, other)
            m.bond(hydrogen, oxygen)

        self.assertEqual(
            canonical_encoding(molecule("a", first)),
            canonical_encoding(molecule("b", second)),
        )

    def test_bond_direction_is_irrelevant(self):
        forward = water()
        reversed_bonds = water()
        reversed_bonds.bonds = [
            Bond(id=b.id, source=b.target, target=b.source, order=b.order)
            for b in reversed_bonds.bonds
        ]
        self.assertEqual(
            canonical_encoding(forward), canonical_encoding(reversed_bonds)
        )

    def test_different_molecules_differ(self):
        distinct = [
            water(),
            hydrogen_peroxide(),
            hydroperoxyl(),
            heavy_water(),
            benzene(),
            methane(),
        ]
        hashes = {content_hash(m) for m in distinct}
        self.assertEqual(len(hashes), len(distinct))

    def test_profile_participates_in_identity(self):
        """
        Spec 10.2: the same graph under a different chemistry model is not the
        same claim, so the profile is inside the hash.
        """
        self.assertNotEqual(
            content_hash(water(), chemistry_profile="pouring:organic-covalent-v0"),
            content_hash(water(), chemistry_profile="something:else-v9"),
        )


# --- stereochemistry -------------------------------------------------------


class TestStereo(unittest.TestCase):
    """
    Spec 4.7 called this disqualifying: without it, a representation cannot
    say which of two enantiomers it means.
    """

    def test_enantiomers_share_formula_and_mass(self):
        left, right = glucose(), glucose(mirrored=True)
        self.assertEqual(formula(left), formula(right))
        self.assertEqual(molar_mass(left), molar_mass(right))

    def test_enantiomers_are_distinguishable(self):
        self.assertNotEqual(
            content_hash(glucose()), content_hash(glucose(mirrored=True))
        )

    def test_mirror_of_mirror_is_the_original(self):
        self.assertEqual(
            content_hash(glucose()), content_hash(mirror(mirror(glucose())))
        )

    def test_geometric_isomers_are_distinguishable(self):
        self.assertEqual(formula(butene()), formula(butene(cis=True)))
        self.assertNotEqual(
            content_hash(butene()), content_hash(butene(cis=True))
        )

    def test_glucose_stereocenters_detected(self):
        self.assertEqual(len(potential_stereocenters(glucose())), 4)

    def test_unannotated_center_is_reported(self):
        bare = glucose()
        bare.stereo_centers = []
        self.assertIn("stereo.unspecified", codes(check_stereo(bare)))

    def test_annotated_center_is_not_reported(self):
        self.assertEqual(check_stereo(glucose()), [])

    def test_severity_is_contextual(self):
        """Spec 4.7: the caller decides, because the module cannot know."""
        bare = glucose()
        bare.stereo_centers = []
        self.assertEqual(check_stereo(bare, severity="error")[0].severity, "error")
        self.assertEqual(check_stereo(bare, severity="info")[0].severity, "info")

    def test_symmetric_molecule_has_no_stereocenters(self):
        self.assertEqual(potential_stereocenters(methane()), [])


# --- balancing -------------------------------------------------------------


class TestBalancer(unittest.TestCase):
    def _resolver(self, route):
        registry = {s.id: s for s in route.species}
        return lambda species_id: registry[species_id]

    def test_solves_a_determined_equation(self):
        route = combustion()
        result = balance(route.reactions[0], self._resolver(route))
        self.assertEqual(result.status, "balanced")
        self.assertEqual(
            result.coefficients,
            {"methane": 1, "dioxygen": 2, "carbon-dioxide": 1, "water": 2},
        )

    def test_renders_the_solved_equation(self):
        route = combustion()
        self.assertEqual(
            balanced_equation(route.reactions[0], self._resolver(route)),
            "methane + 2 dioxygen -> carbon dioxide + 2 water",
        )

    def test_reports_underdetermination_rather_than_guessing(self):
        """
        Aspirin's hydrogen row is exactly twice its oxygen row, so the
        nullspace has dimension 2 and atom counting cannot pick the answer.
        Reporting that beats inventing a coefficient set.
        """
        route = aspirin_synthesis()
        result = balance(route.reactions[0], self._resolver(route))
        self.assertEqual(result.status, "underdetermined")
        self.assertEqual(result.freedom, 2)
        self.assertIsNone(result.coefficients)

    def test_underdetermined_equation_still_validates(self):
        found = validate(
            Document(target_profile=SYNTHESIS, body=aspirin_synthesis())
        )
        self.assertTrue(is_valid(found))
        self.assertIn("reaction.underdetermined", codes(found))

    def test_exact_arithmetic(self):
        """
        Minimal-integer normalisation must be exact. Under floating point a
        coefficient of 3 can arrive as 2.9999999999999996 and scale the whole
        answer off it.
        """
        route = combustion()
        result = balance(route.reactions[0], self._resolver(route))
        for value in result.coefficients.values():
            self.assertIsInstance(value, int)


class TestMultiComponent(unittest.TestCase):
    def test_salt_is_not_a_disconnected_molecule(self):
        salt, parts = sodium_chloride()
        route = Route(
            id="salt",
            name="salt",
            target=salt.id,
            starting_materials=[salt.id],
            species=parts + [salt],
            reactions=[],
        )
        found = validate(Document(target_profile=SYNTHESIS, body=route))
        self.assertNotIn("molecule.disconnected", codes(found))
        self.assertTrue(is_valid(found), [str(d) for d in errors(found)])

    def test_charged_ions_validate(self):
        _, parts = sodium_chloride()
        for ion in parts:
            self.assertEqual(structural_errors(ion), [], ion.name)


# --- units -----------------------------------------------------------------


class TestUnits(unittest.TestCase):
    def test_celsius_is_affine(self):
        self.assertEqual(convert(quantity("90", "degC"), "K").value, Decimal("363.15"))
        self.assertEqual(convert(quantity("0", "degC"), "K").value, Decimal("273.15"))

    def test_gauge_pressure_is_offset_from_absolute(self):
        """1 barg is 2.013 bar absolute, not 1."""
        absolute = convert(quantity("1", "barg"), "bar")
        self.assertAlmostEqual(float(absolute.value), 2.01325, places=4)

    def test_dimension_mismatch_is_rejected(self):
        with self.assertRaises(DimensionMismatch):
            convert(quantity("5", "g"), "mL")

    def test_unknown_unit_is_rejected(self):
        with self.assertRaises(UnknownUnit):
            quantity("5", "furlongs")

    def test_fractions_are_not_interchangeable(self):
        """Converting mass fraction to mole fraction needs composition data."""
        with self.assertRaises(DimensionMismatch):
            convert(quantity("0.5", "massFraction"), "moleFraction")

    def test_round_trip(self):
        original = quantity("2.5", "L")
        self.assertEqual(convert(convert(original, "mL"), "L").value, Decimal("2.5"))

    def test_metamorphic_equivalence(self):
        """Equivalent quantities must normalise identically."""
        self.assertEqual(
            convert(quantity("1", "L"), "mL").value,
            convert(quantity("1000", "mL"), "mL").value,
        )
        self.assertEqual(
            convert(quantity("1", "h"), "s").value,
            convert(quantity("60", "min"), "s").value,
        )

    def test_addition_rejects_affine_units(self):
        """20 degC plus 20 degC is not 40 degC — Celsius has an arbitrary zero."""
        with self.assertRaises(DimensionMismatch):
            add(quantity("20", "degC"), quantity("20", "degC"))

    def test_addition_works_for_absolute_units(self):
        self.assertEqual(
            add(quantity("1", "L"), quantity("500", "mL")).value, Decimal("1.5")
        )

    def test_precision_is_preserved(self):
        """Spec 5.1: "0.100" is a different claim from "0.1"."""
        self.assertEqual(str(quantity("0.100", "mol").value), "0.100")

    def test_conditions_dimension_is_checked(self):
        route = aspirin_synthesis()
        route.reactions[0].conditions.temperature = quantity("90", "mL")
        found = validate(Document(target_profile=SYNTHESIS, body=route))
        self.assertIn("unit.dimension", codes(found))


# --- references and routes -------------------------------------------------


class TestReferenceChecking(unittest.TestCase):
    def test_unresolved_reference(self):
        network = approximate_majority(10, 5)
        network.outputs = ["Z"]
        self.assertIn(
            "ref.unresolved",
            codes(validate(Document(target_profile=CRN, body=network))),
        )

    def test_network_needs_outputs(self):
        network = approximate_majority(10, 5)
        network.outputs = []
        self.assertIn(
            "network.noOutputs",
            codes(validate(Document(target_profile=CRN, body=network))),
        )

    def test_duplicate_ids(self):
        network = approximate_majority(10, 5)
        network.reactions[1].id = "r1"
        self.assertIn(
            "id.duplicate",
            codes(validate(Document(target_profile=CRN, body=network))),
        )

    def test_unreachable_target(self):
        route = combustion()
        route.starting_materials = ["methane"]  # no oxygen supplied
        found = validate(Document(target_profile=SYNTHESIS, body=route))
        self.assertIn("route.unreachableTarget", codes(found))

    def test_reachability_is_a_fixpoint(self):
        """
        The target is produced by a reaction whose own input is unreachable.
        A rule that only asked "is it produced?" would wrongly pass this.
        """
        route = combustion()
        route.starting_materials = []
        found = validate(Document(target_profile=SYNTHESIS, body=route))
        self.assertIn("route.unreachableTarget", codes(found))
        self.assertIn("route.unreachableInput", codes(found))

    def test_byproduct_is_reported(self):
        found = validate(Document(target_profile=SYNTHESIS, body=combustion()))
        self.assertIn("route.unconsumedIntermediate", codes(found))


# --- semantics -------------------------------------------------------------


class TestDeterministicSemantics(unittest.TestCase):
    def test_amplifies_the_majority(self):
        run = simulate_deterministic(approximate_majority(10, 5), t_end=4.0)
        self.assertEqual(majority_verdict(run), "X")

    def test_amplifies_the_other_majority(self):
        run = simulate_deterministic(approximate_majority(5, 10), t_end=4.0)
        self.assertEqual(majority_verdict(run), "Y")

    def test_population_is_conserved(self):
        run = simulate_deterministic(approximate_majority(10, 5), t_end=4.0)
        self.assertAlmostEqual(sum(run.final().values()), 15.0, places=4)

    def test_converges_to_unanimity(self):
        final = simulate_deterministic(
            approximate_majority(10, 5), t_end=6.0
        ).final()
        self.assertAlmostEqual(final["X"], 15.0, places=3)
        self.assertAlmostEqual(final["Y"], 0.0, places=3)


class TestStochasticSemantics(unittest.TestCase):
    def test_reproducible_under_seed(self):
        runs = [
            simulate_stochastic(
                approximate_majority(10, 5, semantics="stochastic"),
                t_end=100.0,
                seed=7,
            ).final()
            for _ in range(2)
        ]
        self.assertEqual(runs[0], runs[1])

    def test_population_is_conserved(self):
        run = simulate_stochastic(
            approximate_majority(10, 5, semantics="stochastic"),
            t_end=100.0,
            seed=1,
        )
        self.assertEqual(sum(run.final().values()), 15.0)

    def test_reaches_unanimity(self):
        final = simulate_stochastic(
            approximate_majority(10, 5, semantics="stochastic"),
            t_end=200.0,
            seed=3,
        ).final()
        self.assertIn(0.0, (final["X"], final["Y"]))

    def test_diverges_from_deterministic_at_low_copy_number(self):
        """
        Spec 6.8's claim. Deterministic always answers X for a 6/4 split;
        stochastic sometimes does not. If this ever passes with zero wrong
        answers, the semantics distinction is decorative and should go.
        """
        wrong = sum(
            majority_verdict(
                simulate_stochastic(
                    approximate_majority(6, 4, semantics="stochastic"),
                    t_end=200.0,
                    seed=seed,
                )
            )
            != "X"
            for seed in range(200)
        )
        self.assertGreater(wrong, 0, "stochastic never diverged")
        self.assertLess(wrong, 100, "stochastic should still favour the majority")

    def test_divergence_shrinks_with_population(self):
        def wrong_rate(x0, y0, runs=200):
            return sum(
                majority_verdict(
                    simulate_stochastic(
                        approximate_majority(x0, y0, semantics="stochastic"),
                        t_end=400.0,
                        seed=seed,
                    )
                )
                != "X"
                for seed in range(runs)
            ) / runs

        self.assertGreater(wrong_rate(6, 4), wrong_rate(60, 40))


# --- encoding --------------------------------------------------------------


class TestEncoding(unittest.TestCase):
    def test_document_envelope(self):
        encoded = Document(
            target_profile=CRN, body=approximate_majority(10, 5)
        ).to_json()
        self.assertEqual(encoded["targetProfile"], CRN)
        self.assertEqual(encoded["schemaVersion"], "0.1.0")
        self.assertEqual(encoded["chemistryProfile"], "pouring:organic-covalent-v0")
        self.assertEqual(encoded["body"]["kind"], "reactionNetwork")

    def test_coefficients_are_exact(self):
        """Spec 5.1: rationals as strings, never floats."""
        encoded = approximate_majority(10, 5).to_json()
        self.assertEqual(
            encoded["reactions"][0]["products"][0]["coefficient"],
            {"numerator": "2", "denominator": "1"},
        )

    def test_json_serialisable(self):
        for body in [autocatalysis(), aspirin_synthesis(), combustion()]:
            profile = CRN if body.kind == "reactionNetwork" else SYNTHESIS
            json.dumps(Document(target_profile=profile, body=body).to_json())

    def test_stereo_survives_encoding(self):
        encoded = glucose().to_json()
        self.assertEqual(len(encoded["stereoCenters"]), 4)

    def test_radical_is_encoded(self):
        encoded = hydroperoxyl().to_json()
        self.assertTrue(any(a.get("radical") for a in encoded["atoms"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
