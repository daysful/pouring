"""
Profile-bound validation.

The spec's central claim is that one reaction set can be read two ways, and
that exactly three rules differ between the readings: balance (6.3), cycles
(6.6-6.7), and kinetics (6.8). This module is where that claim is enforced
rather than asserted — `validate` dispatches on `targetProfile` and the
difference in output is the claim, tested in tests/test_pouring.py.

Structure, stereochemistry, balancing, and units live in their own modules.
This one sequences them and decides, per profile, which apply.
"""

from __future__ import annotations

from .balance import balance, imbalance, is_balanced
from .elements import ValenceOracle, default_oracle
from .ir import (
    BALANCE_SUGGESTED,
    BALANCED,
    CRN,
    ERROR,
    INFO,
    OBSERVATIONAL,
    PARTIAL,
    SYNTHESIS,
    WARNING,
    Diagnostic,
    Document,
    Molecule,
    MultiComponentSpecies,
    ReactionNetwork,
    Route,
    is_structural,
)
from .stereo import check_stereo
from .structure import check_structure
from .units import DimensionMismatch, UnknownUnit


def validate(
    document: Document,
    oracle: ValenceOracle | None = None,
    stereo_severity: str = WARNING,
) -> list[Diagnostic]:
    """
    All checks for one document, in dependency order.

    Graph integrity precedes valence, valence precedes balance, and balance
    precedes anything that trusts a formula — a check run over an unsound
    graph produces confident nonsense.
    """
    body = document.body
    profile = document.target_profile
    oracle = oracle or default_oracle()

    if profile not in (SYNTHESIS, CRN):
        return [
            Diagnostic(
                ERROR, "schema.version", f"unknown targetProfile '{profile}'"
            )
        ]

    found: list[Diagnostic] = []
    found += _check_ids(body)
    found += _check_references(body)
    found += _check_species(body, profile, oracle, stereo_severity)
    found += _check_conditions(body)
    found += _check_balance(body, profile)
    found += _check_kinetics(body, profile)
    found += _check_yields(body)

    if isinstance(body, ReactionNetwork):
        found += _check_network(body)
    else:
        found += _check_route(body)

    return found


def _registry(body) -> dict:
    return {species.id: species for species in body.species}


# --- shared checks ---------------------------------------------------------


def _check_ids(body) -> list[Diagnostic]:
    seen: set[str] = set()
    found = []
    for item in list(body.species) + list(body.reactions):
        if item.id in seen:
            found.append(
                Diagnostic(
                    ERROR,
                    "id.duplicate",
                    f"'{item.id}' declared more than once",
                    path=item.id,
                )
            )
        seen.add(item.id)
    return found


def _check_references(body) -> list[Diagnostic]:
    known = set(_registry(body))
    found = []

    for reaction in body.reactions:
        for participant in list(reaction.reactants) + list(reaction.products):
            if participant.species not in known:
                found.append(
                    Diagnostic(
                        ERROR,
                        "ref.unresolved",
                        f"reaction '{reaction.name}' references unknown "
                        f"species '{participant.species}'",
                        path=reaction.id,
                    )
                )

    if isinstance(body, ReactionNetwork):
        for species_id in list(body.inputs) + list(body.outputs):
            if species_id not in known:
                found.append(
                    Diagnostic(
                        ERROR,
                        "ref.unresolved",
                        f"unknown species '{species_id}'",
                        path=body.id,
                    )
                )

    if isinstance(body, Route):
        for species_id in [body.target] + list(body.starting_materials):
            if species_id not in known:
                found.append(
                    Diagnostic(
                        ERROR,
                        "ref.unresolved",
                        f"unknown species '{species_id}'",
                        path=body.id,
                    )
                )

    return found


def _check_species(
    body,
    profile: str,
    oracle: ValenceOracle,
    stereo_severity: str,
) -> list[Diagnostic]:
    """
    Structural species get the full structure layer. Formal species get an
    explanation of what was skipped (spec 4.8).
    """
    found = []
    known = set(_registry(body))

    for species in body.species:
        if isinstance(species, Molecule):
            found += check_structure(species, oracle)
            found += check_stereo(species, severity=stereo_severity)
            continue

        if isinstance(species, MultiComponentSpecies):
            if not species.components:
                found.append(
                    Diagnostic(
                        ERROR,
                        "species.emptyComponents",
                        f"'{species.name}' declares no components",
                        path=species.id,
                    )
                )
            for component in species.components:
                if component.entity not in known:
                    found.append(
                        Diagnostic(
                            ERROR,
                            "ref.unresolved",
                            f"'{species.name}' references unknown component "
                            f"'{component.entity}'",
                            path=species.id,
                        )
                    )
                if component.ratio <= 0:
                    found.append(
                        Diagnostic(
                            ERROR,
                            "species.invalidRatio",
                            f"component ratio must be positive, got "
                            f"{component.ratio}",
                            path=species.id,
                        )
                    )
            continue

        # Formal species.
        if profile == SYNTHESIS and species.binding is None:
            found.append(
                Diagnostic(
                    ERROR,
                    "species.unresolved",
                    f"'{species.name}' has no structure and no binding; a "
                    f"synthesis cannot be balanced or executed against it",
                    path=species.id,
                )
            )
        else:
            found.append(
                Diagnostic(
                    INFO,
                    "species.formal",
                    f"'{species.name}' has no molecular graph; structural "
                    f"checks skipped",
                    path=species.id,
                )
            )

    return found


def _check_conditions(body) -> list[Diagnostic]:
    """Spec 5.2: every quantity carries a dimension, and they must match."""
    expected = {
        "temperature": "temperature",
        "time": "time",
        "pressure": "pressure",
    }
    found = []

    for reaction in body.reactions:
        if reaction.conditions is None:
            continue

        for field_name, wanted in expected.items():
            value = getattr(reaction.conditions, field_name)
            if value is None:
                continue

            try:
                actual = value.dimension
            except UnknownUnit as error:
                found.append(
                    Diagnostic(
                        ERROR, "unit.unknown", str(error), path=reaction.id
                    )
                )
                continue

            if actual != wanted:
                found.append(
                    Diagnostic(
                        ERROR,
                        "unit.dimension",
                        f"'{field_name}' of '{reaction.name}' is {actual} "
                        f"('{value.unit}'), expected {wanted}",
                        path=reaction.id,
                    )
                )

    return found


def _check_balance(body, profile: str) -> list[Diagnostic]:
    """
    Spec 6.3. Required under synthesis; NOT APPLICABLE under CRN, where
    formal networks violate conservation by design — autocatalysis is a
    standard primitive and no atom accounting permits it.
    """
    if profile != SYNTHESIS:
        return []

    registry = _registry(body)
    found = []

    def resolve(species_id: str):
        return registry[species_id]

    for reaction in body.reactions:
        if reaction.equation_status in (PARTIAL, OBSERVATIONAL):
            continue  # spec 6.4: not every recorded equation is meant to balance

        participants = list(reaction.reactants) + list(reaction.products)
        resolved = [registry.get(p.species) for p in participants]

        if not all(s is not None and is_structural(s) for s in resolved):
            continue  # skipped, not failed — reported by _check_species

        if is_balanced(reaction, resolve):
            # Balanced as written is not the same as uniquely determined. When
            # the nullspace has more than one dimension, other coefficient
            # sets also balance — sometimes absurd ones — and only chemistry
            # rules them out. Worth saying, not worth failing over.
            solved = balance(reaction, resolve)
            if solved.status == "underdetermined":
                found.append(
                    Diagnostic(
                        INFO,
                        "reaction.underdetermined",
                        f"'{reaction.name}' balances as written, but atom "
                        f"counting alone does not determine it: "
                        f"{solved.detail}",
                        path=reaction.id,
                    )
                )
            continue

        left, right = imbalance(reaction, resolve)
        severity = (
            WARNING if reaction.equation_status == BALANCE_SUGGESTED else ERROR
        )

        solved = balance(reaction, resolve)
        hint = ""
        if solved.ok and solved.coefficients:
            rendered = ", ".join(
                f"{registry[k].name}={v}" for k, v in solved.coefficients.items()
            )
            hint = f"; balances at {rendered}"
        elif solved.status == "underdetermined":
            hint = f"; {solved.detail}"

        found.append(
            Diagnostic(
                severity,
                "reaction.unbalanced",
                f"'{reaction.name}' does not balance: "
                f"{_format(left)} vs {_format(right)}{hint}",
                path=reaction.id,
            )
        )

    return found


def _format(tally) -> str:
    from .balance import CHARGE_ROW

    if not tally:
        return "nothing"
    parts = []
    for element, count in sorted(tally.items()):
        if element == CHARGE_ROW:
            parts.append(f"charge{count:+}")
        else:
            parts.append(f"{element}{count}")
    return " ".join(parts)


def _check_kinetics(body, profile: str) -> list[Diagnostic]:
    """
    Spec 6.8. Optional under synthesis; required under CRN, where the rate
    constants are the program.
    """
    if profile != CRN:
        return []

    found = []
    for reaction in body.reactions:
        if reaction.kinetics is None:
            found.append(
                Diagnostic(
                    ERROR,
                    "kinetics.missing",
                    f"'{reaction.name}' has no rate law or rate constant; a "
                    f"network without kinetics has no dynamics to compute with",
                    path=reaction.id,
                )
            )
        elif reaction.kinetics.rate_constant <= 0:
            found.append(
                Diagnostic(
                    ERROR,
                    "kinetics.invalidRate",
                    f"'{reaction.name}' has a non-positive rate constant",
                    path=reaction.id,
                )
            )

    if isinstance(body, ReactionNetwork) and body.semantics is None:
        found.append(
            Diagnostic(
                ERROR,
                "kinetics.semanticsUnspecified",
                "network does not state deterministic or stochastic semantics; "
                "they disagree at low copy number",
                path=body.id,
            )
        )

    return found


def _check_yields(body) -> list[Diagnostic]:
    found = []
    known = set(_registry(body))

    for reaction in body.reactions:
        recorded = reaction.reaction_yield
        if recorded is None:
            continue

        if not 0 <= recorded.value <= 1:
            found.append(
                Diagnostic(
                    ERROR,
                    "reaction.yieldRange",
                    f"yield of '{reaction.name}' is {recorded.value}, outside "
                    f"[0, 1]",
                    path=reaction.id,
                )
            )

        for reference in (recorded.product, recorded.basis):
            if reference not in known:
                found.append(
                    Diagnostic(
                        ERROR,
                        "ref.unresolved",
                        f"yield of '{reaction.name}' references unknown "
                        f"species '{reference}'",
                        path=reaction.id,
                    )
                )

    return found


# --- network checks (CRN reading) ------------------------------------------


def _check_network(network: ReactionNetwork) -> list[Diagnostic]:
    """
    Note what is absent: no cycle check. Catalytic cycles, feedback, and
    oscillators are how a network computes (spec 6.7). Flagging them here
    would forbid computation.
    """
    found = []

    if not network.outputs:
        found.append(
            Diagnostic(
                ERROR,
                "network.noOutputs",
                "network designates no output species, so nothing about it "
                "is checkable",
                path=network.id,
            )
        )

    produced = {p.species for r in network.reactions for p in r.products} | set(
        network.inputs
    )

    for output in network.outputs:
        if output not in produced:
            found.append(
                Diagnostic(
                    WARNING,
                    "network.unreachableOutput",
                    f"output '{output}' is produced by no reaction and is "
                    f"not an input",
                    path=output,
                )
            )

    return found


# --- route checks (synthesis reading) --------------------------------------


def _check_route(route: Route) -> list[Diagnostic]:
    found = []

    # Spec 6.6: reachability is a fixpoint from declared starting materials,
    # not "produced by some reaction" — that reaction's own inputs may be
    # unreachable.
    reachable = set(route.starting_materials)
    changed = True
    while changed:
        changed = False
        for reaction in route.reactions:
            if all(p.species in reachable for p in reaction.reactants):
                for p in reaction.products:
                    if p.species not in reachable:
                        reachable.add(p.species)
                        changed = True

    if route.target not in reachable:
        found.append(
            Diagnostic(
                ERROR,
                "route.unreachableTarget",
                f"target '{route.target}' is not reachable from the declared "
                f"starting materials",
                path=route.id,
            )
        )

    for reaction in route.reactions:
        for participant in reaction.reactants:
            if participant.species not in reachable:
                found.append(
                    Diagnostic(
                        ERROR,
                        "route.unreachableInput",
                        f"'{reaction.name}' needs '{participant.species}', "
                        f"which is never produced or supplied",
                        path=reaction.id,
                    )
                )

    consumed = {p.species for r in route.reactions for p in r.reactants}
    produced = {p.species for r in route.reactions for p in r.products}

    for species in route.species:
        if species.id in consumed or species.id in produced:
            continue
        if species.id in route.starting_materials:
            continue
        found.append(
            Diagnostic(
                WARNING,
                "route.unusedSpecies",
                f"'{species.name}' is declared but neither consumed nor "
                f"produced",
                path=species.id,
            )
        )

    for species_id in produced:
        if species_id != route.target and species_id not in consumed:
            found.append(
                Diagnostic(
                    WARNING,
                    "route.unconsumedIntermediate",
                    f"'{species_id}' is a by-product: produced, never "
                    f"consumed, and not the target",
                    path=species_id,
                )
            )

    if _has_cycle(route):
        found.append(
            Diagnostic(
                ERROR,
                "route.cycle",
                "route contains a cycle; a synthesis plan must terminate",
                path=route.id,
            )
        )

    return found


def _has_cycle(route: Route) -> bool:
    """
    Cycle over the reaction dependency graph: an edge from r to s wherever a
    product of r is a reactant of s.
    """
    edges: dict[str, set[str]] = {r.id: set() for r in route.reactions}

    for first in route.reactions:
        produced = {p.species for p in first.products}
        for second in route.reactions:
            if first.id == second.id:
                continue
            if produced & {p.species for p in second.reactants}:
                edges[first.id].add(second.id)

    visiting: set[str] = set()
    done: set[str] = set()

    def walk(node: str) -> bool:
        if node in visiting:
            return True
        if node in done:
            return False
        visiting.add(node)
        for following in edges[node]:
            if walk(following):
                return True
        visiting.discard(node)
        done.add(node)
        return False

    return any(walk(reaction.id) for reaction in route.reactions)


# --- convenience -----------------------------------------------------------


def errors(diagnostics: list[Diagnostic]) -> list[Diagnostic]:
    return [d for d in diagnostics if d.severity == ERROR]


def warnings(diagnostics: list[Diagnostic]) -> list[Diagnostic]:
    return [d for d in diagnostics if d.severity == WARNING]


def is_valid(diagnostics: list[Diagnostic]) -> bool:
    return not errors(diagnostics)
