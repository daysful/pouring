"""
Worked examples, used by the demo and the tests.

Two families:

  * reaction networks read under `pouring:crn-v0` — programs;
  * molecules and synthesis routes read under `pouring:synthesis-v0` — plans.

`approximate_majority` is the flagship for the first. It is small, it genuinely
computes something, and its deterministic and stochastic readings visibly
disagree, which is what makes the semantics distinction load-bearing rather
than decorative.

The molecules are the conformance corpus of spec 13.2: water against its
near-misses, benzene in both Kekulé forms, and a pair of enantiomers.
"""

from __future__ import annotations

from fractions import Fraction

from .build import molecule
from .ir import (
    BALANCED,
    Component,
    Conditions,
    FormalSpecies,
    Kinetics,
    Molecule,
    MultiComponentSpecies,
    Participant,
    Reaction,
    ReactionNetwork,
    Route,
    Yield,
)
from .units import quantity


def _p(species: str, coefficient: int = 1) -> Participant:
    return Participant(species=species, coefficient=Fraction(coefficient))


# --- molecules -------------------------------------------------------------


def water() -> Molecule:
    """H2O. Two hydrogens on one oxygen."""

    def define(m):
        m.atom("O")
        m.fill()

    return molecule("water", define)


def hydrogen_peroxide() -> Molecule:
    """H2O2 — HO-OH. Same elements as water, different molecule."""

    def define(m):
        first = m.atom("O")
        second = m.atom("O")
        m.bond(first, second)
        m.fill()

    return molecule("hydrogen peroxide", define)


def hydroperoxyl() -> Molecule:
    """
    HO2, the radical the README originally described while calling it water.

    Its middle oxygen carries an unpaired electron, which must be *declared*.
    Spec 4.5: an atom short of its valence is an error, never a silent
    radical — so without `radical=1` this molecule does not validate.
    """

    def define(m):
        first = m.atom("O", radical=1)
        second = m.atom("O")
        m.bond(first, second)
        m.fill()

    return molecule("hydroperoxyl", define)


def benzene(alternate: bool = False) -> Molecule:
    """
    C6H6.

    The ring is the operation an expression tree cannot express. `alternate`
    gives the other Kekulé drawing — the same molecule, and it must
    canonicalise identically (spec 4.2).
    """
    orders = [1, 2, 1, 2, 1, 2] if alternate else [2, 1, 2, 1, 2, 1]

    def define(m):
        ring = m.atoms("C", 6)
        m.ring(ring, orders)
        m.fill()

    return molecule("benzene", define)


# Parities for the four stereocentres of open-chain glucose, in chain order.
# These encode the Fischer convention positionally. Verified CIP descriptors
# are an RDKit responsibility; what matters here is that the representation
# distinguishes an enantiomeric pair at all, which before this it could not.
_D_GLUCOSE_PARITIES = (1, -1, 1, 1)


def glucose(mirrored: bool = False) -> Molecule:
    """
    Open-chain glucose, C6H12O6, with four annotated stereocentres.

    D- and L-glucose have identical atoms, bonds, and formula. The
    specification calls the inability to separate them disqualifying, because
    a system that cannot say which one it means cannot specify what to make.
    """
    parities = [
        -parity if mirrored else parity for parity in _D_GLUCOSE_PARITIES
    ]

    def define(m):
        carbons = m.atoms("C", 6)
        m.chain(carbons)
        m.group(carbons[0], "O", 2)  # aldehyde
        for carbon in carbons[1:]:
            m.hydroxyl(carbon)
        m.fill()
        for carbon, parity in zip(carbons[1:5], parities):
            m.stereocenter(carbon, parity)

    return molecule("L-glucose" if mirrored else "D-glucose", define)


def sodium_chloride() -> tuple[MultiComponentSpecies, list[Molecule]]:
    """
    NaCl as sodium and chloride in 1:1 ratio (spec 4.6).

    Not a covalent molecule with a strange bond, and not one disconnected
    graph — both were available and both are wrong.
    """

    def sodium(m):
        m.atom("Na", charge=1)

    def chloride(m):
        m.atom("Cl", charge=-1)

    parts = [
        molecule("sodium cation", sodium),
        molecule("chloride anion", chloride),
    ]

    salt = MultiComponentSpecies(
        id="sodium-chloride",
        name="sodium chloride",
        components=[
            Component(entity="sodium-cation", ratio=Fraction(1)),
            Component(entity="chloride-anion", ratio=Fraction(1)),
        ],
    )
    return salt, parts


# --- reaction networks -----------------------------------------------------


def approximate_majority(
    x0: int,
    y0: int,
    semantics: str = "deterministic",
    rate: float = 1.0,
) -> ReactionNetwork:
    """
    The approximate-majority network (Angluin et al.; Cardelli).

        X + Y -> 2B     disagreement undecides both
        B + X -> 2X     the undecided join a side
        B + Y -> 2Y

    Given x0 != y0 it amplifies the initial majority to unanimity. Two
    properties matter for the spec:

    It is cyclic by construction — X feeds B feeds X. Under a synthesis
    reading that is `route.cycle`, an error. Here the feedback IS the
    computation.

    Its species are formal. There is no molecular graph for X, Y, or B; they
    are variables, and only a physical realisation (DNA strands, typically)
    would bind them to structure.
    """
    species = [
        FormalSpecies(id="X", name="X"),
        FormalSpecies(id="Y", name="Y"),
        FormalSpecies(id="B", name="B (undecided)"),
    ]

    kinetics = Kinetics(rate_constant=rate)

    reactions = [
        Reaction(
            id="r1",
            name="disagree",
            reactants=[_p("X"), _p("Y")],
            products=[_p("B", 2)],
            kinetics=kinetics,
        ),
        Reaction(
            id="r2",
            name="convert to X",
            reactants=[_p("B"), _p("X")],
            products=[_p("X", 2)],
            kinetics=kinetics,
        ),
        Reaction(
            id="r3",
            name="convert to Y",
            reactants=[_p("B"), _p("Y")],
            products=[_p("Y", 2)],
            kinetics=kinetics,
        ),
    ]

    return ReactionNetwork(
        id="am",
        name="approximate majority",
        species=species,
        reactions=reactions,
        inputs={"X": x0, "Y": y0, "B": 0},
        outputs=["X", "Y"],
        semantics=semantics,
    )


def autocatalysis() -> ReactionNetwork:
    """
    `X -> 2X`. Deliberately violates conservation of mass.

    Under `pouring:synthesis-v0` this is unbalanceable nonsense. Under
    `pouring:crn-v0` it is amplification, a standard primitive, and must not
    raise `reaction.unbalanced` (spec 6.3).
    """
    return ReactionNetwork(
        id="auto",
        name="autocatalysis",
        species=[FormalSpecies(id="X", name="X")],
        reactions=[
            Reaction(
                id="r1",
                name="amplify",
                reactants=[_p("X")],
                products=[_p("X", 2)],
                kinetics=Kinetics(rate_constant=1.0),
            )
        ],
        inputs={"X": 1},
        outputs=["X"],
        semantics="stochastic",
    )


def majority_verdict(trajectory) -> str:
    """What the network decided: the output species that won."""
    final = trajectory.final()
    if abs(final["X"] - final["Y"]) < 1e-6:
        return "tie"
    return "X" if final["X"] > final["Y"] else "Y"


# --- synthesis routes ------------------------------------------------------


def salicylic_acid() -> Molecule:
    """C7H6O3 — benzene with a phenol and a carboxylic acid."""

    def define(m):
        ring = m.atoms("C", 6)
        m.ring(ring, [2, 1, 2, 1, 2, 1])
        m.hydroxyl(ring[0])
        carboxyl = m.group(ring[1], "C")
        m.group(carboxyl, "O", 2)
        m.group(carboxyl, "O")
        m.fill()

    return molecule("salicylic acid", define)


def acetic_anhydride() -> Molecule:
    """C4H6O3 — two acetyl groups sharing a bridging oxygen."""

    def define(m):
        left_methyl = m.atom("C")
        left_carbonyl = m.group(left_methyl, "C")
        m.group(left_carbonyl, "O", 2)
        bridge = m.group(left_carbonyl, "O")
        right_carbonyl = m.group(bridge, "C")
        m.group(right_carbonyl, "O", 2)
        m.group(right_carbonyl, "C")
        m.fill()

    return molecule("acetic anhydride", define)


def aspirin(sabotage: bool = False) -> Molecule:
    """
    C9H8O4 — salicylic acid with the phenol acetylated.

    `sabotage` adds a stray hydrogen as a radical, producing a molecule that
    still validates structurally but no longer balances the equation it
    appears in. It exists to prove the balancer fires rather than merely
    existing.
    """

    def define(m):
        ring = m.atoms("C", 6)
        m.ring(ring, [2, 1, 2, 1, 2, 1])
        ester = m.hydroxyl(ring[0])
        carbonyl = m.group(ester, "C")
        m.group(carbonyl, "O", 2)
        m.group(carbonyl, "C")
        carboxyl = m.group(ring[1], "C")
        m.group(carboxyl, "O", 2)
        m.group(carboxyl, "O")
        m.fill()
        if sabotage:
            m.atom("H", radical=1)

    return molecule("aspirin", define)


def acetic_acid() -> Molecule:
    """C2H4O2."""

    def define(m):
        methyl = m.atom("C")
        carboxyl = m.group(methyl, "C")
        m.group(carboxyl, "O", 2)
        m.group(carboxyl, "O")
        m.fill()

    return molecule("acetic acid", define)


def aspirin_synthesis(sabotage: bool = False) -> Route:
    """
    salicylic acid + acetic anhydride -> aspirin + acetic acid

    C7H6O3 + C4H6O3 -> C9H8O4 + C2H4O2, which balances.
    """
    species = [
        salicylic_acid(),
        acetic_anhydride(),
        aspirin(sabotage=sabotage),
        acetic_acid(),
    ]

    reaction = Reaction(
        id="acetylation",
        name="acetylation of salicylic acid",
        reactants=[_p("salicylic-acid"), _p("acetic-anhydride")],
        products=[_p("aspirin"), _p("acetic-acid")],
        equation_status=BALANCED,
        conditions=Conditions(
            temperature=quantity("90", "degC"),
            time=quantity("20", "min"),
            catalyst="sulfuric acid",
        ),
        reaction_yield=Yield(
            product="aspirin",
            basis="salicylic-acid",
            value=__import__("decimal").Decimal("0.85"),
            kind="isolated",
        ),
    )

    return Route(
        id="aspirin-route",
        name="aspirin from salicylic acid",
        target="aspirin",
        starting_materials=["salicylic-acid", "acetic-anhydride"],
        species=species,
        reactions=[reaction],
    )


def butene(cis: bool = False) -> Molecule:
    """
    2-butene, C4H8 — CH3-CH=CH-CH3, with E/Z configuration declared.

    Geometric isomers, the other half of the v0 stereo model. Like the
    enantiomers, cis and trans share a formula and a mass and must not share a
    canonical form.
    """

    def define(m):
        carbons = m.atoms("C", 4)
        bonds = m.chain(carbons, [1, 2, 1])
        m.fill()
        m.double_bond_stereo(
            bonds[1],
            carbons[0],
            carbons[3],
            "cis" if cis else "trans",
        )

    return molecule("cis-2-butene" if cis else "trans-2-butene", define)


def heavy_water() -> Molecule:
    """
    D2O — water with both hydrogens as deuterium.

    An isotopologue: same graph, same charges, different mass. It must not
    share an identity with water.
    """

    def define(m):
        oxygen = m.atom("O")
        m.bond(oxygen, m.atom("H", isotope=2))
        m.bond(oxygen, m.atom("H", isotope=2))

    return molecule("heavy water", define)


def methane() -> Molecule:
    """CH4."""

    def define(m):
        m.atom("C")
        m.fill()

    return molecule("methane", define)


def dioxygen() -> Molecule:
    """O2."""

    def define(m):
        first = m.atom("O")
        second = m.atom("O")
        m.bond(first, second, 2)

    return molecule("dioxygen", define)


def carbon_dioxide() -> Molecule:
    """CO2 — O=C=O."""

    def define(m):
        carbon = m.atom("C")
        m.group(carbon, "O", 2)
        m.group(carbon, "O", 2)

    return molecule("carbon dioxide", define)


def combustion() -> Route:
    """
    CH4 + 2 O2 -> CO2 + 2 H2O.

    Unlike the aspirin synthesis, this one is uniquely determined: its atom
    matrix has full rank, so the nullspace is one-dimensional and the
    coefficients follow from counting alone. It is the contrast case that
    shows `balance` solving rather than merely checking.
    """
    species = [methane(), dioxygen(), carbon_dioxide(), water()]

    reaction = Reaction(
        id="combustion",
        name="combustion of methane",
        reactants=[_p("methane"), _p("dioxygen", 2)],
        products=[_p("carbon-dioxide"), _p("water", 2)],
        equation_status=BALANCED,
    )

    return Route(
        id="combustion-route",
        name="methane combustion",
        target="carbon-dioxide",
        starting_materials=["methane", "dioxygen"],
        species=species,
        reactions=[reaction],
    )


def approximate_majority_as_route() -> Route:
    """
    The approximate-majority reaction set, read as a synthesis plan.

    Same species, same reactions — a different question asked of them. This is
    the cross-target case from spec 13.2: the differences it produces must be
    exactly the ones the spec predicts, and no others.
    """
    network = approximate_majority(x0=10, y0=5)
    return Route(
        id="am-as-route",
        name="approximate majority, read as a synthesis",
        target="X",
        starting_materials=["X", "Y"],
        species=network.species,
        reactions=network.reactions,
    )
