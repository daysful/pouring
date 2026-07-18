"""
`pouring` — a chemical compiler.

Reference implementation of the IR specified in
about/notes/ir-specification.md.

The CRN target is complete end to end: define a reaction network, validate it,
simulate it under both semantics, check what it computed.

The synthesis target covers structure, stereochemistry, balancing, units, and
route reachability. Valence verdicts come from a `ValenceOracle`, which prefers
RDKit and falls back to a labelled profile lint — chemistry authority is
delegated by design, not reimplemented here.
"""

from .balance import (
    BalanceResult,
    balance,
    balanced_equation,
    is_balanced,
    species_counts,
    suggested_coefficients,
)
from .build import MoleculeBuilder, molecule
from .canonical import (
    canonical_encoding,
    canonical_order,
    content_hash,
    find_rings,
    perceive_aromatic_bonds,
    same_structure,
)
from .elements import (
    MASSES,
    ProfileValenceOracle,
    RDKitValenceOracle,
    ValenceOracle,
    default_oracle,
)
from .ir import (
    BALANCE_SUGGESTED,
    BALANCED,
    CHEMISTRY_V0,
    CRN,
    DETERMINISTIC,
    ERROR,
    INFO,
    NORMALIZE_V0,
    OBSERVATIONAL,
    PARTIAL,
    STOCHASTIC,
    SYNTHESIS,
    WARNING,
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
    is_structural,
)
from .networks import (
    acetic_acid,
    acetic_anhydride,
    approximate_majority,
    approximate_majority_as_route,
    aspirin,
    aspirin_synthesis,
    autocatalysis,
    benzene,
    butene,
    carbon_dioxide,
    combustion,
    dioxygen,
    glucose,
    heavy_water,
    hydrogen_peroxide,
    hydroperoxyl,
    majority_verdict,
    methane,
    salicylic_acid,
    sodium_chloride,
    water,
)
from .simulate import (
    Trajectory,
    simulate,
    simulate_deterministic,
    simulate_stochastic,
)
from .stereo import (
    annotated_centers,
    check_stereo,
    mirror,
    potential_stereocenters,
    stereo_summary,
)
from .structure import (
    check_structure,
    connected_components,
    describe,
    formula,
    hill_formula,
    molar_mass,
    ring_count,
    used_valence,
)
from .units import (
    UNITS,
    DimensionMismatch,
    Quantity,
    UnknownUnit,
    add,
    convert,
    dimension_of,
    quantity,
)
from .validate import errors, is_valid, validate, warnings

__all__ = [
    # profiles and severities
    "CRN",
    "SYNTHESIS",
    "CHEMISTRY_V0",
    "NORMALIZE_V0",
    "DETERMINISTIC",
    "STOCHASTIC",
    "ERROR",
    "WARNING",
    "INFO",
    "BALANCED",
    "BALANCE_SUGGESTED",
    "PARTIAL",
    "OBSERVATIONAL",
    # ir
    "Atom",
    "Bond",
    "Component",
    "Conditions",
    "Diagnostic",
    "Document",
    "DoubleBondStereo",
    "FormalSpecies",
    "Kinetics",
    "Molecule",
    "MultiComponentSpecies",
    "Participant",
    "Reaction",
    "ReactionNetwork",
    "Route",
    "TetrahedralCenter",
    "Yield",
    "is_structural",
    # building
    "MoleculeBuilder",
    "molecule",
    # elements
    "MASSES",
    "ValenceOracle",
    "ProfileValenceOracle",
    "RDKitValenceOracle",
    "default_oracle",
    # structure
    "check_structure",
    "connected_components",
    "describe",
    "formula",
    "hill_formula",
    "molar_mass",
    "ring_count",
    "used_valence",
    # stereo
    "annotated_centers",
    "check_stereo",
    "mirror",
    "potential_stereocenters",
    "stereo_summary",
    # canonical
    "canonical_encoding",
    "canonical_order",
    "content_hash",
    "find_rings",
    "perceive_aromatic_bonds",
    "same_structure",
    # balance
    "BalanceResult",
    "balance",
    "balanced_equation",
    "is_balanced",
    "species_counts",
    "suggested_coefficients",
    # units
    "UNITS",
    "Quantity",
    "UnknownUnit",
    "DimensionMismatch",
    "add",
    "convert",
    "dimension_of",
    "quantity",
    # simulation
    "Trajectory",
    "simulate",
    "simulate_deterministic",
    "simulate_stochastic",
    # validation
    "validate",
    "errors",
    "warnings",
    "is_valid",
    # examples
    "acetic_acid",
    "acetic_anhydride",
    "approximate_majority",
    "approximate_majority_as_route",
    "aspirin",
    "aspirin_synthesis",
    "autocatalysis",
    "benzene",
    "butene",
    "carbon_dioxide",
    "combustion",
    "dioxygen",
    "glucose",
    "heavy_water",
    "methane",
    "hydrogen_peroxide",
    "hydroperoxyl",
    "majority_verdict",
    "salicylic_acid",
    "sodium_chloride",
    "water",
]
