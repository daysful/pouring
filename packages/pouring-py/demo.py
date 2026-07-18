"""
End-to-end demonstration: define, validate, identify, balance, simulate.

    python3 demo.py

No dependencies, no install step.
"""

from pouring import (
    CRN,
    SYNTHESIS,
    Atom,
    Document,
    approximate_majority,
    approximate_majority_as_route,
    aspirin,
    aspirin_synthesis,
    autocatalysis,
    balance,
    balanced_equation,
    benzene,
    butene,
    check_stereo,
    check_structure,
    combustion,
    content_hash,
    convert,
    default_oracle,
    describe,
    formula,
    glucose,
    heavy_water,
    hydrogen_peroxide,
    hydroperoxyl,
    majority_verdict,
    molar_mass,
    potential_stereocenters,
    quantity,
    salicylic_acid,
    simulate_deterministic,
    simulate_stochastic,
    validate,
    water,
)
from pouring.units import DimensionMismatch, add


def rule(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


def report(diagnostics) -> None:
    if not diagnostics:
        print("  (no diagnostics)")
    for diagnostic in diagnostics:
        print(f"  {diagnostic}")


def main() -> None:
    oracle = default_oracle()

    # ---------------------------------------------------------------- 1
    rule("1. The same reaction set, read under both targets")

    network = approximate_majority(x0=10, y0=5, semantics="deterministic")

    print("\napproximate majority:")
    print("    X + Y -> 2B      B + X -> 2X      B + Y -> 2Y")

    print("\n-- as pouring:crn-v0 (a program)")
    report(validate(Document(target_profile=CRN, body=network)))

    print("\n-- as pouring:synthesis-v0 (a plan), identical reactions")
    report(
        validate(
            Document(
                target_profile=SYNTHESIS, body=approximate_majority_as_route()
            )
        )
    )

    print("\n-- and the third difference: strip the rate constants")

    stripped = approximate_majority(x0=10, y0=5)
    for reaction in stripped.reactions:
        reaction.kinetics = None
    stripped.semantics = None

    print("   under crn-v0, where the rate constants are the program:")
    report(
        [
            d
            for d in validate(Document(target_profile=CRN, body=stripped))
            if d.code.startswith("kinetics")
        ]
    )
    print("   under synthesis-v0, where they are optional detail:")
    print("  (no diagnostics)")

    # ---------------------------------------------------------------- 2
    rule("2. Structure: molecules are graphs")

    print()
    for mol in [
        water(),
        hydrogen_peroxide(),
        hydroperoxyl(),
        benzene(),
        glucose(),
        salicylic_acid(),
        aspirin(),
    ]:
        print(f"  {describe(mol)}")

    print(
        "\n  Benzene has a ring because bonds are declared independently of\n"
        "  atom order. A nested expression tree cannot represent one at all."
    )

    print("\n-- valence is checked, and radicals must be declared")
    print(f"   authority: {oracle.name}")

    undeclared = hydroperoxyl()
    undeclared.atoms = [
        Atom(id=a.id, element=a.element, charge=a.charge, radical=0)
        for a in undeclared.atoms
    ]
    print("\n   HO2 with its unpaired electron declared:")
    report([d for d in check_structure(hydroperoxyl(), oracle) if d.severity != "info"])
    print("   the same graph without the declaration:")
    report(check_structure(undeclared, oracle))

    # ---------------------------------------------------------------- 3
    rule("3. Identity: what counts as the same molecule")

    left, right = benzene(), benzene(alternate=True)
    print("\n  benzene drawn two ways (Kekule forms)")
    print(f"    orders   {[b.order for b in left.bonds][:6]}")
    print(f"    orders   {[b.order for b in right.bonds][:6]}")
    print(f"    same     {content_hash(left) == content_hash(right)}")

    d_glucose, l_glucose = glucose(), glucose(mirrored=True)
    print("\n  D-glucose vs L-glucose")
    print(f"    formula  {formula(d_glucose)} / {formula(l_glucose)}")
    print(f"    mass     {molar_mass(d_glucose):.3f} / {molar_mass(l_glucose):.3f}")
    print(f"    centres  {len(potential_stereocenters(d_glucose))}")
    print(f"    same     {content_hash(d_glucose) == content_hash(l_glucose)}")

    print("\n  cis- vs trans-2-butene")
    print(f"    same     {content_hash(butene()) == content_hash(butene(cis=True))}")

    print("\n  water vs heavy water")
    print(f"    mass     {molar_mass(water()):.3f} / {molar_mass(heavy_water()):.3f}")
    print(f"    same     {content_hash(water()) == content_hash(heavy_water())}")

    print(
        "\n  The spec called the inability to separate enantiomers\n"
        "  disqualifying: identical atoms, identical bonds, identical\n"
        "  formula. Now they have different identities."
    )

    print("\n-- an unannotated stereocentre is reported, not assumed away")
    bare = glucose()
    bare.stereo_centers = []
    report(check_stereo(bare)[:2])
    print("  ... and 2 more")

    # ---------------------------------------------------------------- 4
    rule("4. Conservation is profile-bound, not universal")

    print("\nX -> 2X   (autocatalysis: one molecule becomes two)")
    print("\n-- as pouring:crn-v0")
    report(validate(Document(target_profile=CRN, body=autocatalysis())))
    print("  Amplification. No atoms to conserve, so balance is skipped.")

    print("\n-- a synthesis, where balance applies and is solved")
    route = combustion()
    registry = {s.id: s for s in route.species}
    resolve = registry.__getitem__
    print(f"  {balanced_equation(route.reactions[0], resolve)}")
    print(f"  coefficients derived, not assumed: "
          f"{balance(route.reactions[0], resolve).coefficients}")

    print("\n-- and where counting alone is not enough")
    aspirin_route = aspirin_synthesis()
    aspirin_registry = {s.id: s for s in aspirin_route.species}
    result = balance(aspirin_route.reactions[0], aspirin_registry.__getitem__)
    print(f"  aspirin synthesis: {result.status}, freedom {result.freedom}")
    print(
        "  Its hydrogen row is exactly twice its oxygen row, so 0:10:2:11\n"
        "  also balances - and consumes no salicylic acid at all. The\n"
        "  1:1:1:1 a chemist writes is right, but not because of arithmetic."
    )

    print("\n-- the checker fires when it should")
    report(
        [
            d
            for d in validate(
                Document(
                    target_profile=SYNTHESIS,
                    body=aspirin_synthesis(sabotage=True),
                )
            )
            if d.code == "reaction.unbalanced"
        ]
    )

    # ---------------------------------------------------------------- 5
    rule("5. Units carry dimensions")

    print()
    print(f"  90 degC              -> {convert(quantity('90', 'degC'), 'K')}")
    print(f"  1 barg               -> {convert(quantity('1', 'barg'), 'bar')} absolute")
    print(f"  1 L + 500 mL         -> {add(quantity('1', 'L'), quantity('500', 'mL'))}")
    print(f"  0.100 mol            -> {quantity('0.100', 'mol')} (precision kept)")

    for description, thunk in [
        ("5 g -> mL", lambda: convert(quantity("5", "g"), "mL")),
        (
            "massFraction -> moleFraction",
            lambda: convert(quantity("0.5", "massFraction"), "moleFraction"),
        ),
        (
            "20 degC + 20 degC",
            lambda: add(quantity("20", "degC"), quantity("20", "degC")),
        ),
    ]:
        try:
            thunk()
            print(f"  {description:20} -> allowed (unexpected)")
        except DimensionMismatch as error:
            print(f"  {description:20} -> {str(error).replace('[pouring] ', '')}")

    # ---------------------------------------------------------------- 6
    rule("6. Deterministic semantics: mass-action ODEs")

    trajectory = simulate_deterministic(network, t_end=4.0)
    final = trajectory.final()

    print("\n  start   X=10.0  Y=5.0  B=0.0")
    print(f"  end     X={final['X']:.4f}  Y={final['Y']:.4f}  B={final['B']:.4f}")
    print(f"  verdict {majority_verdict(trajectory)}")
    print(f"  total   {sum(final.values()):.4f}  (conserved: started at 15)")

    print("\n  X  ", end="")
    series = trajectory.series("X")
    for index in range(0, len(series), max(1, len(series) // 12)):
        print(f"{series[index]:6.2f}", end="")
    print()

    # ---------------------------------------------------------------- 7
    rule("7. Stochastic semantics: Gillespie over discrete counts")

    print("\n  Same network, same initial majority, 400 independent runs.")
    print("  The deterministic reading always answers X. The stochastic")
    print("  reading is a distribution, and it narrows as counts grow.\n")

    print(f"  {'population':>12}  {'margin':>8}  {'correct':>8}  {'wrong':>7}")
    print(f"  {'-' * 12}  {'-' * 8}  {'-' * 8}  {'-' * 7}")

    for x0, y0 in [(6, 4), (12, 8), (30, 20), (60, 40)]:
        runs = 400
        correct = sum(
            majority_verdict(
                simulate_stochastic(
                    approximate_majority(x0=x0, y0=y0, semantics="stochastic"),
                    t_end=200.0,
                    seed=seed,
                )
            )
            == "X"
            for seed in range(runs)
        )
        rate = 100.0 * correct / runs
        print(
            f"  {x0 + y0:>12}  {x0 - y0:>8}  {rate:>7.1f}%  {100.0 - rate:>6.1f}%"
        )

    print(
        "\n  At 10 molecules the network gets the majority wrong a real\n"
        "  fraction of the time; by 100 it is nearly certain. Two readings\n"
        "  of one network, disagreeing exactly where the spec says they do."
    )


if __name__ == "__main__":
    main()
