"""
Gaut et al. 2026, "A Chemically Defined Synthetic Cell Capable Of Growth And
Replication" — modelled as reaction networks.

The first time `pouring` is pointed at data it did not choose. Everything the
package has been tested against so far was an example selected to demonstrate
the package; these numbers were published by someone else for other reasons,
and either the models reproduce them or they do not.

Two models:

  1. **Selection and competition.** Two synthetic-cell populations differing
     only in the promoter driving alpha-hemolysin, competing for feeder
     liposomes over five generations. Cell + food -> 2 cells is autocatalysis:
     cyclic, non-conserving, kinetics-driven — the CRN target exactly.

  2. **Genome inheritance.** Seven plasmids, no cytoskeleton, no segregation
     machinery, so partitioning at division is random. The paper reports 30%
     of cells holding the complete genome after five generations.

Observed values are read from the figures and main text as published; where a
figure reports a range across replicate experiments, the midpoint is used.
Sources are cited inline against each constant.

    python3 studies/adamala_2026.py
"""

from __future__ import annotations

import math
import random
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pouring import (  # noqa: E402
    CRN,
    Document,
    FormalSpecies,
    Kinetics,
    Participant,
    Reaction,
    ReactionNetwork,
    is_valid,
    simulate_deterministic,
    validate,
)

# --- what the paper reports ------------------------------------------------

# Figure 2c / 4c: fraction of liposomes showing lumen mixing after fusion.
FUSION_SLOW = 0.28  # T7 promoter
FUSION_FAST = 0.45  # T7Max promoter

GENERATIONS = 5

# Figure 4n, Illumina sequencing of T7 vs T7Max abundance after five
# generations of complete cell cycle.
SELECTION_OBSERVED = [
    # (starting fraction of the fast type, observed fraction after 5 gens)
    (0.50, 0.61),
    (0.10, 0.38),
]

# Figure 5f: percentage-point lead of T7Max cells over T7 cells after five
# generations, as feeder liposomes are restricted. Two reciprocal
# fluorophore assignments per condition; midpoint taken.
RESOURCE_OBSERVED = [
    # (feeder level relative to 1x, observed lead in percentage points)
    (1.00, 19.0),   # reported 18 and 20
    (0.50, 25.5),   # reported 25 and 26
    (0.25, 31.0),   # reported 30 and 32
    (0.10, 41.5),   # reported 40 and 43
]

# Figure 3o and 3p.
PLASMID_COUNT = 7
COMPLETE_GENOME_OBSERVED = 0.30       # cells holding all seven plasmids
PER_PLASMID_DETECTION = (0.45, 0.70)  # range across the seven plasmids


def participant(species: str, coefficient: int = 1) -> Participant:
    return Participant(species=species, coefficient=Fraction(coefficient))


# --- model 1: selection and competition ------------------------------------


def competition_network(
    slow: float,
    fast: float,
    feeder: float,
    rate_slow: float,
    rate_fast: float,
) -> ReactionNetwork:
    """
    Two populations feeding from one pool.

        slow + feeder -> 2 slow
        fast + feeder -> 2 fast

    Feeding and division are collapsed into one step because the paper's
    generations are defined that way: a 12 h incubation with feeders, then
    extrusion. What distinguishes the populations is the rate constant, which
    is the whole of the selection claim — same reactions, different kinetics,
    different outcome.
    """
    return ReactionNetwork(
        id="competition",
        name="synthetic cell competition",
        species=[
            FormalSpecies(id="slow", name="T7 cells"),
            FormalSpecies(id="fast", name="T7Max cells"),
            FormalSpecies(id="feeder", name="feeder liposomes"),
        ],
        reactions=[
            Reaction(
                id="feed_slow",
                name="T7 cell feeds and divides",
                reactants=[participant("slow"), participant("feeder")],
                products=[participant("slow", 2)],
                kinetics=Kinetics(rate_constant=rate_slow),
            ),
            Reaction(
                id="feed_fast",
                name="T7Max cell feeds and divides",
                reactants=[participant("fast"), participant("feeder")],
                products=[participant("fast", 2)],
                kinetics=Kinetics(rate_constant=rate_fast),
            ),
        ],
        inputs={"slow": slow, "fast": fast, "feeder": feeder},
        outputs=["slow", "fast"],
        semantics="deterministic",
    )


def run_generations(
    fast_fraction: float,
    feeder_per_generation: float,
    rate_slow: float,
    rate_fast: float,
    generations: int = GENERATIONS,
    duration: float = 1.0,
) -> float:
    """
    Five rounds of: add fresh feeders, incubate, divide, keep a sample.

    Populations are renormalised to a constant total between generations,
    which is what the experiment measures — the paper reports population
    *fractions*, not absolute counts.
    """
    slow = 1.0 - fast_fraction
    fast = fast_fraction

    for _ in range(generations):
        network = competition_network(
            slow, fast, feeder_per_generation, rate_slow, rate_fast
        )
        final = simulate_deterministic(network, t_end=duration, dt=0.002).final()

        slow, fast = final["slow"], final["fast"]
        total = slow + fast
        slow, fast = slow / total, fast / total

    return fast


def fit_rate_ratio(
    start: float,
    observed: float,
    feeder: float = 4.0,
) -> float:
    """
    The rate ratio that reproduces one observation, by bisection.

    Reported separately per experiment on purpose: if the two experiments
    imply different ratios, constant relative fitness is the wrong model and
    averaging them would hide that.
    """
    low, high = 1.0, 20.0
    for _ in range(60):
        middle = (low + high) / 2
        predicted = run_generations(start, feeder, 1.0, middle)
        if predicted < observed:
            low = middle
        else:
            high = middle
    return (low + high) / 2


# --- model 1b: discrete generations with a division cap --------------------


def discrete_generations(
    fast_fraction: float,
    feeder_per_cell: float,
    ratio: float,
    generations: int = GENERATIONS,
) -> float:
    """
    The same competition, but with the experiment's actual structure.

    A cell either captures a feeder in a given round or it does not, and
    division is mechanical extrusion — so a cell divides at most *once* per
    generation however much it ate. Capture probability saturates:

        p_i = 1 - exp(-k_i * feeder_per_cell)

    That cap is the whole difference. When food is plentiful every cell feeds
    and doubles, so the faster type gains nothing; when food is scarce only
    the better competitors feed, and the advantage is expressed in full.
    Scarcity therefore *amplifies* selection, which is what the paper
    reports and what the mass-action network cannot produce.
    """
    slow = 1.0 - fast_fraction
    fast = fast_fraction

    for _ in range(generations):
        probability_slow = 1.0 - math.exp(-feeder_per_cell)
        probability_fast = 1.0 - math.exp(-ratio * feeder_per_cell)

        slow *= 1.0 + probability_slow
        fast *= 1.0 + probability_fast

        total = slow + fast
        slow, fast = slow / total, fast / total

    return fast


def fit_discrete(observations) -> tuple[float, float, float]:
    """
    Grid search over (feeders per cell at 1x, rate ratio).

    Two parameters against four conditions, so it is a real fit rather than
    an interpolation — the model can fail, and the mass-action one did.
    """
    best = (None, None, float("inf"))

    for step in range(1, 121):
        base = step * 0.05
        for ratio_step in range(1, 121):
            ratio = 1.0 + ratio_step * 0.05
            error = 0.0
            for level, observed_lead in observations:
                fast = discrete_generations(0.50, base * level, ratio)
                lead = (fast - (1 - fast)) * 100
                error += (lead - observed_lead) ** 2
            if error < best[2]:
                best = (base, ratio, error)

    return best


# --- model 2: genome inheritance -------------------------------------------


def independent_partition_survival(copies: int, generations: int) -> float:
    """
    Probability a lineage keeps all seven plasmids, if every copy segregates
    independently with even odds.

    Each generation the genome is replicated back to `copies` per plasmid, so
    a plasmid is lost only when all its copies land in the sibling.
    """
    per_plasmid_per_generation = 1.0 - 0.5**copies
    return per_plasmid_per_generation ** (PLASMID_COUNT * generations)


def fit_copy_number(observed: float, generations: int) -> float:
    """Copy number implied by the observed fraction, under independence."""
    per_plasmid = observed ** (1.0 / (PLASMID_COUNT * generations))
    return -math.log2(1.0 - per_plasmid)


def correlated_partition(
    copies: int,
    generations: int,
    bias: float,
    trials: int = 20000,
    seed: int = 0,
) -> tuple[float, float]:
    """
    Monte Carlo partitioning with a shared per-division bias.

    `bias` is how unevenly a division splits lumen contents: 0.5 is a fair
    split, higher means one daughter systematically takes more of everything.
    Because the same split applies to every plasmid, losses become correlated
    — which is the mechanism a bulk-content division would produce, as against
    seven independent coin flips.

    Returns (fraction holding all seven, mean per-plasmid retention).
    """
    rng = random.Random(seed)
    complete = 0
    retained_total = 0

    for _ in range(trials):
        present = [True] * PLASMID_COUNT

        for _ in range(generations):
            # One split fraction for the whole cell, shared by every plasmid.
            share = rng.uniform(1.0 - bias, bias) if bias > 0.5 else 0.5
            for index in range(PLASMID_COUNT):
                if not present[index]:
                    continue
                inherited = sum(
                    1 for _ in range(copies) if rng.random() < share
                )
                if inherited == 0:
                    present[index] = False

        retained_total += sum(present)
        if all(present):
            complete += 1

    return complete / trials, retained_total / (trials * PLASMID_COUNT)


# --- reporting -------------------------------------------------------------


def rule(title: str) -> None:
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


def main() -> None:
    rule("0. The model is a valid pouring document")

    network = competition_network(0.5, 0.5, 4.0, 1.0, 1.6)
    diagnostics = validate(Document(target_profile=CRN, body=network))
    print()
    for diagnostic in diagnostics:
        print(f"  {diagnostic}")
    print(f"\n  valid: {is_valid(diagnostics)}")
    print(
        "\n  Formal species, so no balance check — cell + food -> 2 cells\n"
        "  conserves nothing. Autocatalytic and cyclic, which a synthesis\n"
        "  reading would reject and a network reading requires."
    )

    # ---------------------------------------------------------------- 1
    rule("1. Does one rate ratio explain both selection experiments?")

    print(f"\n  Fusion efficiency reported: {FUSION_SLOW:.0%} (T7) vs "
          f"{FUSION_FAST:.0%} (T7Max)")
    print(f"  Ratio implied by fusion:    {FUSION_FAST / FUSION_SLOW:.2f}\n")

    print(f"  {'start':>7}  {'observed':>9}  {'implied ratio':>14}")
    print(f"  {'-' * 7}  {'-' * 9}  {'-' * 14}")

    fitted = []
    for start, observed in SELECTION_OBSERVED:
        ratio = fit_rate_ratio(start, observed)
        fitted.append(ratio)
        print(f"  {start:>7.0%}  {observed:>9.0%}  {ratio:>14.2f}")

    spread = max(fitted) / min(fitted)
    print(
        f"\n  The two experiments imply ratios differing by {spread:.1f}x.\n"
        f"  Under constant relative fitness one number must explain both,\n"
        f"  so constant relative fitness is not what is happening."
    )

    print("\n  What each fitted ratio predicts for the other experiment:\n")
    print(f"  {'ratio from':>12}  {'predicts 50% ->':>16}  {'predicts 10% ->':>16}")
    print(f"  {'-' * 12}  {'-' * 16}  {'-' * 16}")
    for (start, _), ratio in zip(SELECTION_OBSERVED, fitted):
        first = run_generations(0.50, 4.0, 1.0, ratio)
        second = run_generations(0.10, 4.0, 1.0, ratio)
        print(f"  {start:>11.0%}   {first:>15.0%}   {second:>15.0%}")
    print(f"  {'observed':>12}  {0.61:>15.0%}   {0.38:>15.0%}")

    # ---------------------------------------------------------------- 2
    rule("2. Does restricting food widen the lead, as reported?")

    print(
        "\n  Figure 5f: the T7Max lead grows as feeders are withheld. An\n"
        "  explicit shared food pool should reproduce the direction without\n"
        "  being told to.\n"
    )

    ratio = fitted[0]
    print(f"  {'feeder':>8}  {'observed lead':>14}  {'modelled lead':>14}")
    print(f"  {'-' * 8}  {'-' * 14}  {'-' * 14}")

    modelled = []
    for level, observed_lead in RESOURCE_OBSERVED:
        fast = run_generations(0.50, 4.0 * level, 1.0, ratio)
        lead = (fast - (1 - fast)) * 100
        modelled.append(lead)
        print(f"  {level:>7.2f}x  {observed_lead:>13.1f}  {lead:>13.1f}")

    direction = "yes" if modelled[-1] > modelled[0] else "NO"
    print(f"\n  Lead widens as food shrinks: {direction}")

    if direction == "NO":
        print(
            "\n  The mass-action network gets this backwards, and the reason\n"
            "  is structural rather than a bad parameter. In a reaction\n"
            "  network, less substrate simply means less reaction: fewer\n"
            "  feeding events, so less time for the rate difference to act,\n"
            "  so *less* divergence. No rate constant reverses that."
        )

    # ---------------------------------------------------------------- 2b
    rule("2b. What the experiment's structure actually is")

    print(
        "\n  Division is mechanical extrusion: a cell divides at most once\n"
        "  per generation however much it fed. That cap is not expressible\n"
        "  as a reaction — it is a discrete-time population rule. Adding it,\n"
        "  and letting capture probability saturate:\n"
    )

    base, discrete_ratio, error = fit_discrete(RESOURCE_OBSERVED)

    print(f"  fitted feeders per cell at 1x : {base:.2f}")
    print(f"  fitted rate ratio             : {discrete_ratio:.2f}")
    print(f"  residual (sum of squares)     : {error:.1f}\n")

    print(f"  {'feeder':>8}  {'observed lead':>14}  {'modelled lead':>14}")
    print(f"  {'-' * 8}  {'-' * 14}  {'-' * 14}")

    discrete_leads = []
    for level, observed_lead in RESOURCE_OBSERVED:
        fast = discrete_generations(0.50, base * level, discrete_ratio)
        lead = (fast - (1 - fast)) * 100
        discrete_leads.append(lead)
        print(f"  {level:>7.2f}x  {observed_lead:>13.1f}  {lead:>13.1f}")

    widens = "yes" if discrete_leads[-1] > discrete_leads[0] else "NO"
    print(f"\n  Lead widens as food shrinks: {widens}")
    print(
        "\n  When food is plentiful every cell feeds and doubles, so the\n"
        "  faster type gains nothing. When food is scarce only the better\n"
        "  competitors feed at all. Scarcity does not slow selection down —\n"
        "  it is what lets selection happen."
    )

    # ---------------------------------------------------------------- 3
    rule("3. Is genome inheritance consistent with independent segregation?")

    implied = fit_copy_number(COMPLETE_GENOME_OBSERVED, GENERATIONS)
    print(
        f"\n  Observed: {COMPLETE_GENOME_OBSERVED:.0%} of cells hold all "
        f"{PLASMID_COUNT} plasmids after {GENERATIONS} generations."
    )
    print(f"  Copy number implied under independent segregation: {implied:.1f}")

    rounded = max(1, round(implied))
    predicted_complete = independent_partition_survival(rounded, GENERATIONS)
    per_plasmid = (1 - 0.5**rounded) ** GENERATIONS

    print(f"\n  At {rounded} copies per plasmid, independence predicts:")
    print(f"    complete genome     {predicted_complete:>6.1%}   "
          f"(observed {COMPLETE_GENOME_OBSERVED:.0%})")
    print(f"    per-plasmid retention {per_plasmid:>6.1%}   "
          f"(observed {PER_PLASMID_DETECTION[0]:.0%}-"
          f"{PER_PLASMID_DETECTION[1]:.0%})")

    midpoint = sum(PER_PLASMID_DETECTION) / 2
    independent_all_seven = midpoint**PLASMID_COUNT
    print(
        f"\n  The consistency check that matters: if the seven plasmids "
        f"segregated\n  independently, then holding all seven would occur at "
        f"the product of\n  the individual rates."
    )
    print(f"    mean per-plasmid detection   {midpoint:.0%}")
    print(f"    product over {PLASMID_COUNT} plasmids        "
          f"{independent_all_seven:.1%}")
    print(f"    actually observed            {COMPLETE_GENOME_OBSERVED:.0%}")
    print(
        f"\n  {COMPLETE_GENOME_OBSERVED / independent_all_seven:.0f}x more "
        f"complete genomes than independence allows. Inheritance is\n"
        f"  strongly correlated: cells tend to keep most plasmids or lose "
        f"many,\n  rather than losing them one at a time."
    )

    print(
        "\n  Searching copy number and split bias for a pair that reproduces\n"
        "  BOTH reported numbers — 30% complete and ~57% per plasmid. Either\n"
        "  alone is easy to hit; together they constrain the mechanism.\n"
    )

    print(f"  {'copies':>7}  {'bias':>6}  {'complete':>9}  {'per-plasmid':>12}  "
          f"{'error':>7}")
    print(f"  {'-' * 7}  {'-' * 6}  {'-' * 9}  {'-' * 12}  {'-' * 7}")

    best = None
    for copies in (3, 4, 5, 6, 8):
        for bias in (0.50, 0.75, 0.90, 0.98):
            complete, retention = correlated_partition(
                copies, GENERATIONS, bias, trials=6000
            )
            error = (complete - COMPLETE_GENOME_OBSERVED) ** 2 + (
                retention - midpoint
            ) ** 2
            if best is None or error < best[0]:
                best = (error, copies, bias, complete, retention)

    for copies in (3, 4, 5, 6, 8):
        for bias in (0.50, 0.98):
            complete, retention = correlated_partition(
                copies, GENERATIONS, bias, trials=6000
            )
            error = (complete - COMPLETE_GENOME_OBSERVED) ** 2 + (
                retention - midpoint
            ) ** 2
            marker = " <-" if best and (copies, bias) == best[1:3] else ""
            print(
                f"  {copies:>7}  {bias:>6.2f}  {complete:>8.1%}  "
                f"{retention:>11.1%}  {error:>7.3f}{marker}"
            )

    if best:
        _, copies, bias, complete, retention = best
        print(
            f"\n  Closest: {copies} copies with a {bias:.2f} split bias — "
            f"{complete:.0%} complete,\n  {retention:.0%} per plasmid, against "
            f"{COMPLETE_GENOME_OBSERVED:.0%} and {midpoint:.0%} observed."
        )
        print(
            "\n  Neither extreme fits. Independent segregation cannot produce\n"
            "  30% complete genomes at these per-plasmid rates, and a fair\n"
            "  even split cannot produce the losses. What fits is uneven bulk\n"
            "  partitioning: daughters that inherit most of the lumen keep\n"
            "  everything, and daughters that inherit little lose several\n"
            "  plasmids at once."
        )

    print(
        "\n  Correlated loss is what a cell dividing its bulk contents does.\n"
        "  Independent loss is what seven separately segregated plasmids\n"
        "  would do — and the paper says there is no segregation machinery."
    )

    # ---------------------------------------------------------------- 4
    rule("4. What the models found")

    print("""
  ON THE PAPER

  1. The two selection experiments are not consistent with each other
     under constant relative fitness. Fitting the 1:1 result predicts
     15% for the 9:1 experiment, which reported 38%; fitting the 9:1
     result predicts 83% for the 1:1 experiment, which reported 61%.
     Something frequency-dependent is happening, or one measurement is
     further off than its error bars suggest.

  2. Reproductive advantage is far smaller than fusion advantage. Fusion
     efficiency differs by 1.6x; the per-generation reproductive ratio
     needed to explain the outcomes is 1.06-1.23. Feeding better does not
     translate proportionally into leaving more offspring.

  3. Genome inheritance is not independent, by a wide margin. At the
     reported per-plasmid rates, independent segregation allows about 2%
     of cells to hold all seven plasmids. The paper reports 30% — some
     14x more. Losses must be strongly correlated, which points to uneven
     bulk partitioning of lumen contents rather than plasmid-by-plasmid
     segregation. That is consistent with the paper's own statement that
     there is no segregation machinery, and it is not a claim the paper
     makes.

  ON POURING

  4. The population model is an ordinary CRN document and validates as
     one: formal species, no balance check, cycles required rather than
     forbidden. The profile-bound rules were bound correctly for a case
     nobody designed them around.

  5. But the mass-action network gets the resource result backwards, and
     no rate constant fixes it. The experiment caps division at once per
     generation regardless of feeding, and a cap is not a reaction. This
     is a real limit of the formalism as implemented, not a modelling
     slip: `simulate.py` only integrates mass action, though `Kinetics`
     already declares `rateLaw` as an open field.

  CAVEATS

  Observed values are read from published figures, not underlying data.
  The discrete model reproduces the direction of the resource effect and
  the rough magnitude, but not its shape — it peaks at 0.25x where the
  reported lead is still climbing at 0.1x. Two parameters against four
  points is a weak fit and should be treated as indicative.
""")


if __name__ == "__main__":
    main()
