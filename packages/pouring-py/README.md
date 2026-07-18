# pouring — Python reference implementation

Reference implementation of the [IR specification](../../about/notes/ir-specification.md).

Per the [implementation-language decision](../../about/notes/implementation-language.md),
Python is where claims about chemistry live, and this package owns validation
and semantics. TypeScript is *allocated* the eDSL, schema tooling, editor, and
browser playground — none of which exist yet. Until they do, this is the only
implementation.

```bash
python3 demo.py                                  # end to end, ~1s
python3 -m unittest discover -s tests -v         # 99 tests, ~0.8s
```

No dependencies, no install step.



## Modules

| Module | Covers |
|---|---|
| `ir.py` | the IR as data, with canonical JSON encoding |
| `codec.py` | reading documents back — decode, round trip, diagnostics |
| `smiles.py` | SMILES import, with Kekulisation |
| `build.py` | molecular graph construction — atoms, bonds, rings, elaboration |
| `elements.py` | atomic and isotope masses; the valence oracle seam |
| `structure.py` | formula, mass, connectivity, rings, the valence rule |
| `stereo.py` | tetrahedral parity, E/Z, stereocentre detection |
| `canonical.py` | Morgan labelling, Kekulé normalisation, content hashing |
| `balance.py` | exact rational nullspace balancing |
| `units.py` | dimensions, affine and gauge conversions |
| `validate.py` | profile-bound checks — the same set read two ways |
| `simulate.py` | mass-action ODE (RK4) and stochastic (Gillespie) |
| `networks.py` | the conformance corpus |



## Where chemistry authority lives

`elements.py` separates two kinds of data, because they deserve different
treatment. **Atomic masses are facts** — published numbers, safe to hard-code.
**Permitted valences are a model**, and the specification concedes its own rule
is not precise enough for two implementations to agree.

So valence goes behind a `ValenceOracle` protocol. `default_oracle()` prefers
RDKit and falls back to a labelled profile lint; diagnostics name which
authority answered. Callers depend on the protocol, never the table. Install
the real thing with `pip install -e '.[chemistry]'`.



## Getting structures in and out

```python
from pouring import from_smiles, dumps, loads

aspirin = from_smiles('CC(=O)Oc1ccccc1C(=O)O')   # C9H8O4, 1 ring
```

Before SMILES import the only way to state a molecule was to write
graph-construction code, which is fine for a corpus and hopeless for a person.
Parsed structures canonicalise identically to hand-built ones — two independent
construction paths converging on one identity, which is the cross-check that
neither is quietly wrong.

The encoding also reads back, which until recently it did not: `to_json`
existed and nothing parsed it, making an interchange format that did not
interchange. `dumps(round_trip(d)) == dumps(d)` now holds for every body type,
stereochemistry included — a configuration lost in transit would make the
molecule that arrives a different claim from the one that was sent.

Malformed input produces diagnostics rather than exceptions. A bad document is
ordinary input to a compiler, and a caller loading someone else's file deserves
a list of what is wrong with it, not a traceback at the first bad field.

## The claim under test

The specification asserts one reaction set can be read as either a synthesis
plan or a chemical program, and that **exactly three rules differ**: balance,
cycles, and kinetics.

`TestCrossTarget::test_no_fourth_difference` enforces that rather than
asserting it — it runs the approximate-majority network through both profiles
and fails if they diverge in any unlicensed way. If it goes red, the shared-IR
claim is wrong.

```
-- as pouring:crn-v0 (a program)
  info    species.formal [X]: 'X' has no molecular graph; structural checks skipped

-- as pouring:synthesis-v0 (a plan), identical reactions
  error   species.unresolved [X]: 'X' has no structure and no binding ...
  error   route.cycle [am-as-route]: route contains a cycle; a synthesis plan must terminate
```



## Three findings worth knowing

**Identity is not the graph.** Benzene's two Kekulé drawings are the same
molecule, so rings whose bond orders alternate are encoded with a symbolic
aromatic bond and hash identically. Meanwhile D- and L-glucose share a formula,
a mass, and every bond, and must *not* hash identically. Both now hold.

**Balance is not always determined.** The aspirin synthesis balances as written
and its nullspace still has dimension 2, because its hydrogen row is exactly
twice its oxygen row. `0 C7H6O3 + 10 C4H6O3 → 2 C9H8O4 + 11 C2H4O2` also
balances while consuming no salicylic acid. The balancer reports this rather
than picking one, which is what spec 6.3 means by "the answer is a choice".

**Semantics genuinely disagree.** Approximate majority under the deterministic
reading always returns the initial majority. Under the stochastic reading it
does not:

```
  population    margin   correct    wrong
          10         2     89.2%    10.8%
          20         4     94.5%     5.5%
          50        10     99.8%     0.2%
```

Low copy number is where molecular computation operates, which is why the IR
requires a network to declare which reading it means.



## Not implemented

SMILES export (import only); chirality markers `@`/`@@` and directional bonds
in SMILES input, which are *reported* rather than silently dropped; molecule
templates and ports; atom-mapped reaction transforms; resonance,
tautomer, and protonation normalisation (so a hash means "same annotated
graph", not "same chemical entity"); full CIP priority, axial chirality, and
non-tetrahedral coordination; the procedure and execution layers; and every
backend — XDL, Autoprotocol, DNA strand displacement.

See specification §14.
