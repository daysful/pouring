# Studies

Applications of `pouring` to data it did not choose.

Everything in `tests/` and `networks.py` is an example selected to demonstrate
the package. These are not. The point of a study is that it can fail, and one
of them does.

```bash
python3 studies/adamala_2026.py
```



## adamala_2026 — synthetic cell growth, selection, and inheritance

Models Gaut et al. 2026, *A Chemically Defined Synthetic Cell Capable Of Growth
And Replication* (Gaut, Deich, Cash, Hoog, Engelhart, Adamala; U. Minnesota).

Liposome-encapsulated synthetic cells carrying a 90 kbp genome across seven
plasmids, feeding by fusion with "feeder" liposomes, dividing by extrusion, and
undergoing five generations of selection between two populations that differ
only in the promoter driving alpha-hemolysin.

Stripped of biology, the system is:

```
cell + feeder -> 2 cell
```

Autocatalysis. Cyclic, non-conserving, and driven entirely by the rate
difference between the two populations — the CRN target's three profile-bound
rules, all bound correctly for a case they were not designed around.

### Findings about the paper

**The two selection experiments do not agree with each other.** Fitting a
constant relative fitness to the 1:1 result predicts 15% for the 9:1
experiment, which reported 38%. Fitting the 9:1 result predicts 83% for the
1:1 experiment, which reported 61%. One rate ratio cannot produce both.

**Reproductive advantage is much smaller than fusion advantage.** Fusion
efficiency differs 1.6x (45% vs 28%); the per-generation reproductive ratio
needed to explain the outcomes is 1.06–1.23. Feeding better does not translate
proportionally into more offspring.

**Genome inheritance is strongly correlated, not independent.** At the reported
per-plasmid detection rates (~57% mean), independent segregation permits about
2% of cells to hold all seven plasmids. The paper reports **30%** — roughly 14x
more. The losses must be correlated, which points to uneven bulk partitioning
of lumen contents rather than plasmid-by-plasmid segregation. A search over
copy number and split bias lands near 8 copies with a strongly uneven split.

This is consistent with the paper's own statement that the cells have no
cytoskeleton and no DNA segregation mechanism, and it is not a claim the paper
makes.

### Finding about pouring

**The mass-action network gets the resource result backwards.** The paper
reports the faster population's lead *widening* as feeders are withheld — 19 to
41.5 percentage points from 1x down to 0.1x. The CRN predicts the opposite,
22 down to 3.5, and no rate constant reverses it: in a reaction network, less
substrate simply means less reaction, so less divergence.

The experiment's structure is the missing piece. Division is mechanical
extrusion, so a cell divides at most **once per generation** however much it
fed. When food is plentiful every cell feeds and doubles and the faster type
gains nothing; when food is scarce only the better competitors feed at all.
Scarcity does not slow selection — it is what lets selection happen.

A cap of that kind is not a reaction. It is a discrete-time population rule,
and it sits outside the formalism as implemented: `simulate.py` integrates mass
action only, though `Kinetics.rate_law` is already an open field.

That is the useful part of a negative result — it names a missing feature
rather than a wrong number.

### Caveats

Observed values are read from published figures, not underlying data; where a
figure reports a range across reciprocal replicate experiments, the midpoint is
used. The discrete model reproduces the direction and rough magnitude of the
resource effect but not its shape, peaking at 0.25x where the reported lead is
still climbing at 0.1x. Two parameters against four points is a weak fit and is
indicative only.
