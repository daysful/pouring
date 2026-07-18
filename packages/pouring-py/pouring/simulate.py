"""
Semantics for reaction networks (spec 6.8).

A network with kinetics admits two standard readings, and they disagree:

    deterministic   mass-action ODEs over continuous concentrations
    stochastic      a CTMC over discrete molecule counts, via Gillespie

The disagreement is not a rounding artefact. At low copy number the stochastic
reading can settle on a different answer than the deterministic one, and low
copy number is exactly where molecular computation operates. That is why the
IR requires a network to declare which reading it means.

Pure standard library on purpose: this must run with no install step. RK4 is
ample for small non-stiff mass-action systems. Stiff or large networks want
scipy's solvers, which is where the Python-side dependency starts earning its
keep.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from .ir import ReactionNetwork


@dataclass
class Trajectory:
    species: list[str]
    times: list[float] = field(default_factory=list)
    states: list[list[float]] = field(default_factory=list)
    semantics: str = ""

    def record(self, time: float, state: list[float]) -> None:
        self.times.append(time)
        self.states.append(list(state))

    def final(self) -> dict[str, float]:
        return dict(zip(self.species, self.states[-1]))

    def series(self, species_id: str) -> list[float]:
        index = self.species.index(species_id)
        return [state[index] for state in self.states]

    def dominant(self) -> str:
        """The species with the largest final amount."""
        final = self.final()
        return max(final, key=final.get)


def _initial_state(network: ReactionNetwork, order: list[str]) -> list[float]:
    return [float(network.inputs.get(species, 0.0)) for species in order]


# --- deterministic ---------------------------------------------------------


def _derivatives(
    network: ReactionNetwork,
    order: list[str],
    state: list[float],
) -> list[float]:
    index = {species: i for i, species in enumerate(order)}
    delta = [0.0] * len(order)

    for reaction in network.reactions:
        rate = reaction.kinetics.rate_constant
        for participant in reaction.reactants:
            concentration = state[index[participant.species]]
            rate *= concentration ** float(participant.coefficient)

        for species, change in reaction.net_change().items():
            delta[index[species]] += float(change) * rate

    return delta


def simulate_deterministic(
    network: ReactionNetwork,
    t_end: float,
    dt: float = 0.001,
    sample_every: int = 50,
) -> Trajectory:
    """Mass-action ODEs, integrated with fixed-step RK4."""
    order = network.species_ids()
    trajectory = Trajectory(species=order, semantics="deterministic")

    state = _initial_state(network, order)
    time = 0.0
    trajectory.record(time, state)

    steps = int(t_end / dt)
    for step in range(steps):
        k1 = _derivatives(network, order, state)
        mid1 = [s + dt / 2 * d for s, d in zip(state, k1)]
        k2 = _derivatives(network, order, mid1)
        mid2 = [s + dt / 2 * d for s, d in zip(state, k2)]
        k3 = _derivatives(network, order, mid2)
        end = [s + dt * d for s, d in zip(state, k3)]
        k4 = _derivatives(network, order, end)

        state = [
            max(0.0, s + dt / 6 * (a + 2 * b + 2 * c + d))
            for s, a, b, c, d in zip(state, k1, k2, k3, k4)
        ]
        time += dt

        if (step + 1) % sample_every == 0:
            trajectory.record(time, state)

    if trajectory.times[-1] < time:
        trajectory.record(time, state)

    return trajectory


# --- stochastic ------------------------------------------------------------


def _propensity(reaction, counts: dict[str, int]) -> float:
    """
    Mass-action propensity over discrete counts: k times the number of
    distinct reactant combinations available.

        A + B -> ...     k * n_A * n_B
        2A    -> ...     k * n_A * (n_A - 1) / 2
    """
    value = reaction.kinetics.rate_constant
    for participant in reaction.reactants:
        available = counts.get(participant.species, 0)
        needed = int(participant.coefficient)
        if available < needed:
            return 0.0
        value *= math.comb(available, needed)
    return value


def simulate_stochastic(
    network: ReactionNetwork,
    t_end: float,
    seed: int | None = None,
    max_events: int = 1_000_000,
) -> Trajectory:
    """Gillespie's direct method over discrete molecule counts."""
    rng = random.Random(seed)
    order = network.species_ids()
    trajectory = Trajectory(species=order, semantics="stochastic")

    counts = {
        species: int(round(network.inputs.get(species, 0))) for species in order
    }
    time = 0.0
    trajectory.record(time, [float(counts[s]) for s in order])

    for _ in range(max_events):
        propensities = [_propensity(r, counts) for r in network.reactions]
        total = sum(propensities)

        if total <= 0.0:
            break  # no reaction can fire; the network has settled

        time += -math.log(rng.random()) / total
        if time > t_end:
            break

        threshold = rng.random() * total
        cumulative = 0.0
        chosen = network.reactions[-1]
        for reaction, propensity in zip(network.reactions, propensities):
            cumulative += propensity
            if cumulative >= threshold:
                chosen = reaction
                break

        for species, change in chosen.net_change().items():
            counts[species] = counts.get(species, 0) + int(change)

        trajectory.record(time, [float(counts[s]) for s in order])

    trajectory.record(t_end, [float(counts[s]) for s in order])
    return trajectory


# --- dispatch --------------------------------------------------------------


def simulate(network: ReactionNetwork, t_end: float, **kwargs) -> Trajectory:
    """Runs the semantics the network declares."""
    if network.semantics == "stochastic":
        return simulate_stochastic(network, t_end, **kwargs)
    return simulate_deterministic(network, t_end, **kwargs)
