# pouring IR — specification

**Status:** draft, v0. Nothing here is stable yet.
**Scope:** this document specifies the *chemical* IR only — identity, materials,
and transformation intent. The procedure, platform, and execution layers are
declared in §7–§8 but deliberately not yet specified.

> **Revision note.** This draft corrects four defects in the previous one: an
> unreachable `route.orphan` diagnostic, a missing pressure dimension, an
> incorrect route-reachability rule, and an overbroad claim about yield
> composition. It also narrows the valence rule from a general theory to a
> profile-bound lint, and separates the entities that were previously blurred
> together under "compound".
>
> It further admits a **second target**: computing *with* molecules, not only
> making them (§3.1). Three rules previously stated as universal turn out to
> hold only for synthesis — balance (§6.3), acyclicity (§6.6), and kinetics
> being out of scope (§6.8) — and species are no longer required to have a
> molecular structure at all (§4.8).

---

## 0. Why an IR at all

`pouring` is a language whose output is eventually consumed by machines that
move real matter. Between the text a human writes and the instructions a
platform executes there needs to be one well-specified artefact both sides
agree on.

That artefact is this IR. It exists so that **many frontends** can produce it —
a TypeScript or Python eDSL, `.pouring` surface syntax, a SMILES importer, a
structure editor — **many backends** can consume it, and **checking happens
once**, in the middle, rather than being reimplemented badly in every tool.

The IR is data, not code, with a canonical JSON encoding (§10). This document
is normative; implementations are not.

---

## 1. Design principles

**1.1 — Layers are never conflated.** Six questions, six answers (§2). The
original design used one primitive, `react`, for two of them.

**1.2 — Structure is a netlist, not a tree.** Molecules are graphs. Bonds are
declared independently of atom order, so a ring is a bond that closes back. A
nested expression tree cannot represent benzene; a netlist can.

**1.3 — The IR is explicit; the source may be implicit.** Chemists write
implicit hydrogens. The IR stores every atom. Frontends saturate during
elaboration.

**1.4 — Conservation laws are free checks, where they apply.** Atom and charge
balance are objective and cheap, requiring no chemical knowledge beyond
counting. They hold under the synthesis target and are *not applicable* under
the CRN target, where formal networks violate conservation by design (§6.3).

**1.5 — Scope is declared, never assumed.** Every document names the chemistry
profile it was written against (§3). Anything outside that profile is a
`chemistry.unsupported` error, never a silently wrong answer.

**1.6 — Necessary is not sufficient.** Every check here establishes a necessary
condition for a route to work. None establishes sufficiency, and no accumulation
of them ever will (§12).

---

## 2. Entity ontology

The previous draft used *compound*, *fragment*, *species*, and *pour* without
distinguishing them. Those terms name genuinely different things, and conflating
them is a reliable source of error.

| Entity | Meaning |
|---|---|
| `Element` | A chemical element type — oxygen |
| `Atom` | One node in a molecular graph — *this* oxygen, with these neighbours |
| `Molecule` | A connected molecular graph with charge and stereochemistry |
| `MultiComponentSpecies` | Salts, solvates, hydrates, co-crystals — components with ratios |
| `Species` | Chemical identity under a stated protonation/tautomer model |
| `Substance` | Bulk composition, possibly impure or multicomponent |
| `MaterialBatch` | An actual physical batch: amount, purity, lot, phase, provenance |
| `Aliquot` | An identified portion of a batch |
| `MoleculeTemplate` | A graph with attachment sites — not itself a molecule |
| `ReactionEquation` | A balanced statement of what converts into what |
| `Transfer` | An *operation* that moves material |

Two consequences worth stating plainly:

**An element is not an atom.** `element('O')` names oxygen the element; it
cannot identify a particular oxygen in a molecule. Two carbons in one molecule
share an element but differ in identity, neighbours, and stereochemical role.
Frontends should say `atom({ element: 'O' })`.

**`Pour` was two things at once.** An aliquot is a material object; pouring is
an operation, and only one of several ways to transfer. The IR entity is
`Aliquot` (§5.4); `Transfer` belongs to the procedure layer (§7). The authoring
verb `pour` may remain — it is the project's name — but it elaborates to an
aliquot, not to an operation.

---

## 3. Document envelope

Every IR document carries:

```json
{
  "schemaVersion": "0.1.0",
  "targetProfile": "pouring:synthesis-v0",
  "chemistryProfile": "pouring:organic-covalent-v0",
  "normalizationProfile": "pouring:normalize-v0"
}
```

`schemaVersion` governs the encoding. `targetProfile` governs what the document
*means* (§3.1). `chemistryProfile` governs what chemistry is representable and
how valence, aromaticity, and stereochemistry are interpreted. All participate
in identity hashing (§10.2) — the same graph under a different profile is not
necessarily the same claim.

### 3.1 Two targets over one core

`pouring` serves two purposes that share a representation but not their rules:

| | `pouring:synthesis-v0` | `pouring:crn-v0` |
|---|---|---|
| Purpose | make a molecule | compute with molecules |
| Species | must be structural | may be formal (§4.8) |
| Mass/charge balance | required | **not applicable** |
| Kinetics | optional | **required** |
| Cycles | error in a route | expected and essential |
| Semantics | a plan to execute | a dynamical system to evolve |
| Backends | XDL, Autoprotocol | DNA strand displacement, simulation |

**What is shared is the reaction set**: species, reactions, coefficients,
conditions, diagnostics, and — where species are structural — the whole
structure layer. What differs is the *interpretation* laid over it. A synthesis
route and a chemical program are two readings of the same underlying object,
which is the entire argument for one IR rather than two projects.

Three rules below would be silently wrong if applied across both, and are
therefore profile-bound rather than universal: balance (§6.3), cycles
(§6.6–6.7), and kinetics (§6.8).

### 3.2 The v0 chemistry profile

`pouring:organic-covalent-v0` covers:

- a restricted element set;
- formal charges and isotopes;
- covalent bonds of order 1, 2, 3;
- closed-shell molecules, plus simple monoradicals (§4.5);
- tetrahedral and double-bond stereochemistry *once §4.7 lands*.

It explicitly does **not** cover: coordinate/dative bonds, zero-order and
quadruple bonds, multicentre and haptic organometallic bonding, ionic
association within a single graph, delocalisation not reducible to a supported
Kekulé form, carbenes and other multi-electron open-shell species, or
macromolecules.

Encountering any of these must raise `chemistry.unsupported`. Producing a
plausible-looking valence result for chemistry the profile does not model is
worse than refusing.

---

## 4. Structure layer

### 4.1 Atom

| Field | Type | Required | Meaning |
|---|---|---|---|
| `id` | string | yes | unique within the document |
| `element` | string | yes | element symbol; must be in the profile |
| `charge` | integer | yes | formal charge; `0` if neutral |
| `radical` | integer | no | unpaired electrons; default `0` (see §4.5) |
| `isotope` | integer | no | mass number, when a specific isotope is meant |
| `label` | string | no | diagnostics only; carries no semantics |

### 4.2 Bond

| Field | Type | Required | Meaning |
|---|---|---|---|
| `id` | string | yes | unique within the document |
| `from` | AtomId | yes | endpoint |
| `to` | AtomId | yes | endpoint; must differ from `from` |
| `order` | 1 \| 2 \| 3 | yes | bond order in whole units |

Bonds are undirected. `from`/`to` carry no chemical meaning, and implementations
must treat `{a,b}` and `{b,a}` as identical — this is a conformance requirement
(§13).

**Aromaticity is not a bond order.** Aromatic systems are stored in Kekulé form
so valence arithmetic stays exact; aromaticity is a *perception pass*, derived
on demand under a named model, never stored.

This is a real choice with a real cost: aromaticity is model-dependent, and
equivalent Kekulé forms of benzene must be recognised as identical during
canonicalisation or they will hash differently (§10.2). Deriving is still
preferable to admitting an "aromatic" bond kind whose qualifying rings are
undefined.

### 4.3 Port

An intentionally open valence, so a fragment can be wired into a larger
structure later.

| Field | Type | Required | Meaning |
|---|---|---|---|
| `id` | string | yes | unique within the document |
| `name` | string | yes | unique within its template |
| `atom` | AtomId | yes | the atom carrying the open valence |
| `width` | 1 \| 2 \| 3 | yes | bond-order units reserved |

**Ports belong to templates, not molecules.** A `MoleculeTemplate` has ports; a
`Molecule` does not. Once elaboration finishes, a manufacturable molecule with
open ports is an error, not a partially-specified compound. Ports will
eventually need more than a width — permitted bond kind, stereochemical
consequence, and whether attachment removes a capping hydrogen or leaving group.

### 4.4 Molecule

`{ kind: "molecule", id, name, atoms[], bonds[] }`, connected.

Disconnected graphs are **not** molecules. Salts, hydrates, solvates, and
co-crystals use `MultiComponentSpecies` (§4.6) rather than being smuggled in as
one disconnected molecule.

### 4.5 The valence rule — a profile lint, not a theory

For every atom *a* within `pouring:organic-covalent-v0`:

```
used(a) = Σ order(b)  over bonds b incident to a
        + Σ width(p)  over ports p on a          [templates only]
        + radical(a)

V(a)    = valences permitted for element(a) under the profile,
          adjusted for charge(a)
```

Then exactly one of: `used(a) ∈ V(a)` → satisfied; `used(a) > max(V(a))` →
`valence.exceeded`; otherwise → `valence.floating`.

An atom short of its valence is an error, not a silent radical. Unpaired
electrons must be *declared*, never inferred from an omission.

**Known limits of this rule — it is a useful lint for a restricted organic
subset and nothing more:**

- Radical electrons do not universally behave as bond-order units. This works
  for simple monoradicals (methyl: 3 bonds + 1 radical = 4) and fails for
  carbenes, where singlet and triplet states differ in ways one integer cannot
  express. Such species are outside the profile.
- Permitted valence depends on element, charge, bonding context, aromaticity,
  and documented special cases. Hypervalent and organometallic species need
  richer semantics than a lookup table.
- "Adjusted for charge" is **not specified precisely enough for two independent
  implementations to agree.** Until it is, the profile binds to a named
  reference implementation (§13.2), and the conformance corpus — not this prose
  — is what actually pins the behaviour.

### 4.6 Multi-component species

```
{ kind: "multiComponentSpecies", id, name,
  components: [ { entity: EntityRef, ratio: Rational } ] }
```

Sodium chloride is sodium and chloride in 1:1 ratio, not a covalent molecule
with a strange bond. Hydrates and solvates carry their solvent as a component
with an explicit ratio.

### 4.7 Stereochemistry — unimplemented and disqualifying

**This remains the largest gap and it blocks any claim of exact chemical
identity.** D-glucose and L-glucose have identical atoms, bonds, and formula,
as do the other fourteen aldohexoses. A representation that cannot separate
them cannot specify what to manufacture.

Planned, in rough priority order:

1. tetrahedral chirality (ordered neighbour list plus parity);
2. E/Z double-bond configuration;
3. explicitly *unknown* and explicitly *racemic* stereochemistry — different
   claims, both distinct from "unannotated";
4. enhanced stereo groups (relative vs absolute configuration);
5. axial chirality — allenes, cumulenes, atropisomerism;
6. square-planar, trigonal-bipyramidal, octahedral coordination.

Items 5–6 are outside the v0 profile and should raise
`chemistry.unsupported` rather than be silently dropped.

`stereo.unspecified` needs **contextual severity**, not a fixed one. For an
exact API it is an error; for a declared racemate it is intentional; for an
intermediate where stereochemistry is irrelevant it is acceptable; for an
imported historical record it is simply unknown. Severity is therefore set by
the consuming profile, not by this document.

### 4.8 Formal species

Not every species has a structure, and requiring one was an assumption inherited
from the synthesis target.

```
{ kind: "formalSpecies", id, name, binding?: SpeciesRef }
```

A formal species is an identity with no molecular graph. In `pouring:crn-v0`
this is the normal case: the reaction network `X + Y → Z` is a *program*, and
`X`, `Y`, `Z` are its variables. They acquire structure only when the network is
compiled to a physical implementation — DNA strands, typically — at which point
`binding` records what realised them.

This is the same relationship a register has to a physical wire, and it is why
the abstract network is worth specifying separately from any realisation of it.

Formal species are useful in `pouring:synthesis-v0` too, for reagents referenced
but not drawn. There they must resolve before a document is considered complete
(`species.unresolved`), because a synthesis whose reagents have no structure
cannot be balanced, costed, or executed.

Structural and formal species are interchangeable wherever a `SpeciesRef`
appears. Checks that require structure — balance, formula, molar mass — are
skipped with `species.formal` at info severity rather than failing.

---

## 5. Material layer

### 5.1 Numerics

**JSON `number` is not the canonical numeric type.** Binary floating point
breaks equality, hashing, canonical serialisation, exact stoichiometry, and
significant figures. The IR encodes:

- `DecimalValue` — a string, `"0.100"`, preserving precision as written;
- `Rational` — `{ numerator, denominator }` as strings, for exact coefficients;
- `IntegerValue` — a string.

Stoichiometric coefficients are `Rational`. Measured quantities are
`DecimalValue`.

### 5.2 Dimensions and units

| Dimension | Units |
|---|---|
| amount | `mol`, `mmol` |
| mass | `g`, `mg`, `kg` |
| volume | `L`, `mL` |
| temperature | `K`, `degC` |
| time | `s`, `min`, `h` |
| concentration | `M` |
| **pressure** | `Pa`, `kPa`, `bar`, `atm`, `mmHg` |
| flow rate | `mL/min`, `L/h` |
| stirring rate | `rpm` |
| length | `m`, `cm`, `mm` |
| density | `g/mL` |
| fraction | `massFraction`, `moleFraction`, `volumeFraction` |
| equivalents | `equiv` |
| wavelength | `nm` |
| potential | `V` |
| current | `A` |
| acidity | `pH` |

Operations between mismatched dimensions are rejected (`unit.dimension`).

Two conversions are **not** simple multiplication and must be handled
specially: Celsius is affine, not scalar; and gauge pressure differs from
absolute and must be distinguished rather than assumed.

### 5.3 Quantity

```
{ value?: DecimalValue,
  range?: { min?: DecimalValue, max?: DecimalValue },
  unit: UnitRef,
  uncertainty?: DecimalValue,
  qualifier?: "exact" | "approximate" | "lessThan" | "greaterThan" }
```

A measurement is not a number. Charging "about 5 mL", "5.00 mL", and "at least
5 mL" are different claims, and a system that flattens them loses information a
chemist deliberately recorded.

**The unit is mandatory.** `pour(water, 1)` is not a quantity of anything: one
mole and one gram of water are different objects, and unit confusion is a
well-documented way to destroy expensive things.

### 5.4 Aliquot

```
{ kind: "aliquot", id, species: SpeciesRef, amount: Quantity,
  batch?: MaterialBatchRef }
```

**Quantities do not create materials on their own.** Converting between mass,
amount, and volume may require molar mass, concentration, purity, density,
temperature, phase, and composition. An aliquot that names only a species and a
number is under-determined for anything but the simplest case — which is why
`MaterialBatch` exists even though v0 does not yet model it fully.

---

## 6. Transformation layer

Previously called the "process" layer. That was a misnomer: this layer models
*chemical transformation intent*, not physical process. Process is §7.

### 6.1 ReactionEquation

Renamed from `Reaction`, because that is what it actually is — an equation, not
an experiment and not a procedure.

| Field | Type | Required | Meaning |
|---|---|---|---|
| `kind` | `"reactionEquation"` | yes | |
| `id`, `name` | string | yes | |
| `reactants` | Participant[] | yes | |
| `products` | Participant[] | yes | |
| `conditions` | Conditions | no | |
| `equationStatus` | enum | yes | see §6.4 |

`Participant` is `{ species: SpeciesRef, coefficient: Rational }`.

**Equation coefficients and experimental equivalents are different quantities.**
An equation says `2 H₂ + O₂ → 2 H₂O`; an experiment says "1.00 equiv substrate,
1.20 equiv reagent, 0.05 equiv catalyst". The first is exact and stoichiometric;
the second is a charging decision. They do not share a type — experimental use
belongs to `MaterialUse` in the procedure layer (§7).

**Solvent and catalyst are references, not strings.** Free strings cannot
support inventory checking, quantity calculation, hardware execution, hazard
analysis, or provenance. v0 permits strings only as an import shim, flagged
`reaction.unresolvedMaterial`.

### 6.2 Atom mapping

An equation without atom mapping cannot support reaction templates,
retrosynthesis, isotope tracking, mass provenance, or reaction-centre analysis.
Three distinct types are therefore needed:

- `ReactionEquation` — species and coefficients, balanced;
- `AtomMappedReaction` — plus an atom-to-atom correspondence;
- `ReactionTransform` — an atom-mapped graph rewrite usable as a template.

Only the first is specified in v0. The other two are required before any
retrosynthesis work begins.

### 6.3 Balance

**Profile-bound: required under `pouring:synthesis-v0`, not applicable under
`pouring:crn-v0`.**

Summing atom counts across reactants weighted by coefficient must equal the sum
over products, per element; total formal charge must likewise match. Violations
are `reaction.unbalanced` and `reaction.charge`.

Conservation was described in §1.4 as a free check, and under the synthesis
target it is. Under the CRN target it is **wrong**: formal reaction networks
routinely and deliberately violate it. Autocatalysis — `X → X + X` — is a
standard CRN primitive and the basis of amplification, yet no accounting of
atoms permits it. A formal network is a program over abstract species, and
abstract species have no atoms to conserve.

Balance therefore applies exactly when every participant is structural. It is
skipped, not failed, otherwise. When such a network is later compiled to a
physical realisation, the *realisation* must balance even though the network it
implements did not — conservation reappears at the layer where matter does.

This is a nullspace computation over the integer atom-and-charge matrix, `Aν =
0`. **Solving it is not the whole problem.** A balancer must also handle: no
solution; a nullspace of dimension greater than one (the equation is
underdetermined and the answer is a *choice*); positivity and sign constraints
separating reactants from products; minimal-integer normalisation; coefficients
fixed by the author; and optional balancing species (H₂O, H⁺, OH⁻, electrons)
whose admissibility depends on declared conditions.

Exact rational arithmetic is required (§5.1). Floating point makes minimal
integer normalisation unreliable.

### 6.4 Equation status

Not every recorded equation is meant to balance. A literature record may
deliberately omit by-products; an executable mass ledger may not.

```
equationStatus: "authorSpecified" | "balanced" | "balanceSuggested"
              | "partial" | "observational"
```

`reaction.unbalanced` is an error only for `balanced`; for `partial` and
`observational` it is suppressed, and for `balanceSuggested` it is a warning
carrying the suggested coefficients as a fix.

### 6.5 Yield

A single number is insufficient. Yield requires at minimum: the product it
refers to; the basis or limiting reagent; whether it is isolated, assay, or
calculated; the value with uncertainty; the analytical method; and the scale.
Conversion and selectivity are separate quantities and must not be folded in.

```
{ product: SpeciesRef, basis: SpeciesRef,
  kind: "isolated" | "assay" | "calculated",
  value: Quantity, method?: string }
```

Absent yield means unknown, never quantitative.

### 6.6 Route

**Profile-bound: `pouring:synthesis-v0` only.** A route is a *plan* — it has a
target, it terminates, and it is acyclic. The CRN reading of the same reaction
set is §6.7.

```
{ kind: "route", id, name, target: SpeciesRef,
  startingMaterials: SpeciesRef[],
  species[], reactions[] }
```

**Starting materials are explicit.** The previous draft inferred them — "species
produced by no reaction" — which silently reclassified every stray or misspelt
species as procurable, and made `route.orphan` unreachable by construction.
Inference here trades a loud error for a quiet one.

**The network is a hypergraph, not a graph.** A reaction consumes several
species and produces several, which is not a binary edge. Represent it as a
directed hypergraph, or equivalently as a bipartite graph of species nodes and
reaction nodes. A *selected* route is acyclic; the retrosynthesis search space
it is drawn from is an AND/OR structure — a target is made by alternative
reactions (OR), each needing all its precursors (AND).

**Reachability is a fixpoint, not a lookup.** The previous rule — target is
produced by some reaction — is insufficient, since that reaction's inputs may
themselves be unreachable. Correctly:

1. mark declared starting materials reachable;
2. mark a reaction reachable when *all* its reactants are reachable;
3. mark that reaction's products reachable;
4. repeat to fixpoint;
5. the target must be reachable.

**Yields do not generally compose by multiplication.** For a strictly linear
1:1 route the overall fractional yield is the product of step yields. For
convergent, branching, recycled, or split-stream routes it is a material-flow
problem requiring scale propagation and limiting-material selection at each
node. The previous claim that yields simply "compose along" a route was wrong
and is withdrawn.

### 6.7 Reaction network

**Profile-bound: `pouring:crn-v0`.**

```
{ kind: "reactionNetwork", id, name,
  species[], reactions[],
  inputs: [ { species: SpeciesRef, initial: Quantity } ],
  outputs: SpeciesRef[] }
```

The same species-and-reaction set as a route, read as a dynamical system rather
than a plan. The differences are not cosmetic:

**Cycles are the mechanism, not an error.** `route.cycle` must not be raised
here. Catalytic cycles, feedback loops, oscillators, and bistable switches are
how a network computes; forbidding cycles would forbid computation. A network
has no "target" and need not terminate — a chemical oscillator that halts has
failed.

**Inputs and outputs replace starting materials and target.** A network is
evaluated by setting initial concentrations and observing designated output
species over time. This is the testbench, and without it nothing about a network
is checkable.

**Reaction networks are Turing-universal.** Arbitrary CRNs can be implemented by
DNA strand displacement, which is what makes this target physically real rather
than a metaphor — the compilation chain runs abstract network → strand
displacement gates → orderable DNA sequences. Existing work to interoperate with
rather than duplicate: Visual DSD, Nuskell, Peppercorn, and CRN++ (an imperative
language that compiles to reaction networks — the closest prior art to what
`pouring` would be for this target).

### 6.8 Kinetics

**Profile-bound: optional under `pouring:synthesis-v0`, required under
`pouring:crn-v0`.**

```
{ rateLaw: "massAction" | "michaelisMenten" | ...,
  rateConstant: Quantity,
  reversible?: boolean }
```

Under the synthesis target, kinetics are useful but rarely decisive — a chemist
mostly wants to know what forms, not how fast.

Under the CRN target, **the rate constants are the program**. A network's
computation *is* its dynamics; two networks with identical reactions and
different rate constants compute different functions. Declaring kinetics out of
scope, as the previous draft did, would have made the computational target
unrepresentable.

A network with kinetics admits two standard semantics, and the IR must say which
is intended because they disagree:

- **deterministic** — mass-action ODEs over continuous concentrations;
- **stochastic** — a continuous-time Markov chain over discrete molecule counts,
  simulated by the Gillespie algorithm.

The distinction matters at low copy number, where stochastic behaviour diverges
sharply from the deterministic limit, and low copy number is exactly where
molecular computation operates.

**With no hardware, simulation is the whole verification story** — and unlike
the synthesis target, it is a complete one. A chemical program can be run,
observed, and checked against its claimed function entirely in software. That
makes the CRN target the shorter path to something that demonstrably works.

---

## 7. Procedure layer — declared, not specified

**This is the largest gap between what `pouring` currently is and what it claims
to be.** A reaction equation states an intended chemical conversion. It does not
say what a laboratory physically does, and no amount of detail in §6 will make
it do so.

An automation-capable system needs operations: `Transfer`, `Dispense`,
`AddSolid`, `Heat`, `Cool`, `Stir`, `Wait`, `Measure`, `Sample`, `Filter`,
`Wash`, `Separate`, `Extract`, `Evaporate`, `Distill`, `Crystallize`, `Dry`,
`Quench`, `Clean`.

Each needs explicit inputs and outputs, a target vessel, a control objective and
tolerance, temporal semantics, preconditions, completion criteria, failure
modes, and resource requirements.

**The intended strategy is to compile to XDL and Autoprotocol, not to invent a
competing vocabulary.** Both already define hardware-independent chemical
operations, and XDL already compiles to robotic platforms. `pouring`'s
contribution is the typed layer above them, not a replacement for them.

### 7.1 Elaboration-time versus runtime control flow

An embedded DSL creates a semantic hazard that must be resolved here rather than
left to each frontend. In a TypeScript or Python eDSL:

```ts
if (temperature > 50) { cool(); }
```

executes *while constructing the IR*. It does not encode a runtime conditional —
it decides, once, at elaboration, what IR to emit. A runtime condition must be
an explicit IR node:

```ts
when(sensor('reactor.temperature').greaterThan(quantity('50 degC')),
     cool({ vessel: reactor }));
```

Likewise host-language `for` unrolls at elaboration, while `repeatUntil(...)`
must survive into the IR as control flow. Every hardware DSL confronts this;
`pouring` must name which constructs are elaboration-time and which are runtime,
in the specification, not in each frontend's documentation.

---

## 8. Platform binding and execution — roadmap

Two further IRs are anticipated below the procedure layer:

**Execution IR** — a plan bound to one equipment capability graph: bound
vessels, pumps and sensors, scheduled operations, resource locks, interlocks,
recovery actions. Checked against real equipment limits and dry-runnable.

**Runtime trace** — facts rather than intentions: commands issued and accepted,
measurements observed, tolerances exceeded, operator approvals, completions,
failures, material actually transferred.

The full pipeline this document sits at the top of:

```
            eDSLs / .pouring / structure editor
                          ↓ elaboration
                  Chemical IR  ← this document
                          ↓ normalisation, validation
                          |
         ┌────────────────┴────────────────┐
         ↓ synthesis                       ↓ crn
   Procedure IR (§7)              Reaction network (§6.7)
         ↓ scheduling                      ↓ compilation
   Execution IR                    Strand-displacement gates
         ↓                                 ↓
  XDL / Autoprotocol              DNA sequences / simulation
         ↓                                 ↓
   Runtime trace                    Observed trajectory
```

Both branches terminate in a record of what actually happened, and both are
checked against what was intended. The right-hand branch is reachable today
without hardware, because its verification step is simulation; the left-hand
branch is not.

Keeping these separate is what stops the identity model filling up with pumps,
and stops the execution model believing a molecule name is enough to run
hardware.

---

## 9. Diagnostics

```
{ severity: "error" | "warning" | "info",
  code: string,
  message: string,
  path?: JsonPointer,
  sourceSpan?: SourceSpan,
  related?: RelatedLocation[],
  fixes?: SuggestedFix[] }
```

Codes are stable; messages are not. Tools match on `code`.

A bare `subject: string` was too weak for nested documents: a diagnostic must
locate itself with a JSON Pointer *and*, where the frontend supplied one, a span
back into authoring source. Without source mapping, `valence.floating` can name
an atom id but not the line or the drawn object that produced it.

| Code | Severity | Condition |
|---|---|---|
| `schema.version` | error | unknown or unsupported `schemaVersion` |
| `chemistry.unsupported` | error | construct outside the declared profile |
| `id.duplicate` | error | repeated id within a document |
| `ref.unresolved` | error | reference to a non-existent entity |
| `species.formal` | info | structural check skipped — species has no graph (§4.8) |
| `species.unresolved` | error² | formal species with no binding, under a synthesis profile |
| `element.unknown` | error | element not in the profile |
| `molecule.empty` | error | no atoms |
| `molecule.disconnected` | error | use `MultiComponentSpecies` instead |
| `bond.dangling` | error | endpoint references an unknown atom |
| `bond.selfLoop` | error | `from` equals `to` |
| `bond.duplicate` | error | two bonds share an atom pair — raise the order |
| `valence.floating` | error | unsatisfied valence, undeclared (§4.5) |
| `valence.exceeded` | error | bonding exceeds maximum valence |
| `port.dangling` | error | port references an unknown atom |
| `port.duplicate` | error | duplicate port name |
| `port.onMolecule` | error | open port on a non-template molecule |
| `stereo.unspecified` | *profile* | potential stereocentre unannotated (§4.7) |
| `reaction.unbalanced` | error¹ | atom counts differ across the arrow |
| `reaction.charge` | error¹ | total formal charge differs |
| `reaction.unresolvedMaterial` | warning | solvent/catalyst given as a bare string |
| `reaction.openPorts` | error | a reaction species is a template, not a molecule |
| `reaction.yieldRange` | error | yield outside [0,1] |
| `unit.unknown` | error | unrecognised unit |
| `unit.dimension` | error | dimension mismatch |
| `route.unreachableTarget` | error | target not reachable from starting materials (§6.6) |
| `route.unreachableInput` | error | a reaction's reactant is never produced or supplied |
| `route.unusedSpecies` | warning | declared but neither consumed nor produced |
| `route.unconsumedIntermediate` | warning | produced, not the target, never consumed |
| `route.cycle` | error³ | the selected route contains a cycle |
| `network.noOutputs` | error | reaction network designates no output species |
| `network.unreachableOutput` | warning | an output species is produced by no reaction |
| `kinetics.missing` | error⁴ | reaction has no rate law or rate constant |
| `kinetics.semanticsUnspecified` | error⁴ | network does not state deterministic or stochastic |

¹ Balance codes: error under `pouring:synthesis-v0`, not raised under
`pouring:crn-v0` (§6.3). Further suppressed or downgraded per `equationStatus`
(§6.4).
² Raised only under a synthesis profile; formal species are normal under
`pouring:crn-v0` (§4.8).
³ Raised only under a synthesis profile. Cycles are expected in a reaction
network and must not be flagged (§6.7).
⁴ Raised only under `pouring:crn-v0`, where kinetics are required (§6.8).

`route.orphan` from the previous draft is **removed**: under its own definition
of starting materials it could never fire. The four `route.*` codes above replace
it with conditions that can actually occur.

---

## 10. Canonical encoding and identity

### 10.1 Encoding

JSON. Field order is insignificant. Ids are opaque strings, unique within a
document, carrying no meaning — implementations must not parse them.

Water:

```json
{
  "schemaVersion": "0.1.0",
  "chemistryProfile": "pouring:organic-covalent-v0",
  "kind": "molecule",
  "id": "m0",
  "name": "water",
  "atoms": [
    { "id": "a0", "element": "O", "charge": 0 },
    { "id": "a1", "element": "H", "charge": 0 },
    { "id": "a2", "element": "H", "charge": 0 }
  ],
  "bonds": [
    { "id": "b0", "from": "a0", "to": "a1", "order": 1 },
    { "id": "b1", "from": "a0", "to": "a2", "order": 1 }
  ]
}
```

Two hydrogens on one oxygen. The README's original example placed the count on
the wrong atom and described HO₂; under §4.5 that form is a `valence.floating`
error rather than a silent mistake, which is the entire argument for a checked
IR.

### 10.2 Identity is harder than canonical numbering

Deterministic atom ids do **not** yield chemical identity. Before any hash is
meaningful, the normalisation profile must fix policies for: alternate Kekulé
forms, resonance, tautomerism, protonation state, salt disconnection, isotopes,
formal-charge normalisation, stereochemistry, enhanced stereo groups, and
implicit versus explicit hydrogens. InChI itself separates normalisation,
canonicalisation, and serialisation for exactly these reasons.

An identifier must therefore commit to its semantics, not just its graph:

```
hash( schemaVersion, chemistryProfile, normalizationProfile,
      canonicalStructureEncoding )
```

Hashing the visible JSON graph alone would make two documents equal that assert
different things.

Documents should additionally carry established identifiers rather than
replacing them:

```json
"identifiers": [
  { "type": "inchi", "value": "..." },
  { "type": "inchikey", "value": "..." },
  { "type": "canonicalSmiles", "value": "..." },
  { "type": "molfile", "value": "..." }
]
```

Canonicalisation remains deferred to v1. "Morgan/InChI-style" was a gesture at a
solution, not a specification, and is not treated as one here.

---

## 11. Safety

A system that may eventually move matter cannot treat safety as one more
validator pass. Independent boundaries, in order:

1. **Chemical validation** — structure, charge, quantities, references consistent
2. **Procedure validation** — operations have defined inputs, outputs, completion criteria
3. **Capability validation** — equipment can satisfy temperature, pressure, volume, flow, material compatibility
4. **Hazard review** — incompatibilities, runaway risk, toxic release, pressure generation, waste
5. **Dry run** — equipment graph and schedule exercised without material
6. **Human authorisation** — a qualified operator approves an executable plan
7. **Runtime interlocks** — hardware independently enforces limits and safe states
8. **Execution trace** — requested and actual values permanently recorded

This IR addresses layer 1 only. **No output of this specification is
authorisation to execute anything physically**, and tooling must not present a
clean validation as if it were.

---

## 12. Non-goals

**Mechanism is out of scope — but kinetics are not.** The previous draft
grouped them, which was a mistake: they are different claims. Mechanism
*explains* a transformation — electron pushing, transition states — and remains
out of scope. Kinetics *describe observable dynamics*, and are required under
the CRN target where they constitute the program itself (§6.8). `pouring`
describes what transforms into what, how fast, and under which conditions —
never why.

**Feasibility is out of scope, permanently.** A balanced, valence-valid,
well-typed route can still fail in a flask. Every check here is a *necessary*
condition and none is sufficient — not because the checks are immature, but
because chemistry does not compose the way logic gates do. Protecting groups
exist precisely because functional groups interfere non-locally.

"Balanced and valence-valid" is not "will work", and is not "safe". This
limitation belongs in the README, not buried here.

---

## 13. Conformance

### 13.1 Test classes

Prose such as "valence adjusted by charge" cannot make two implementations
agree. The corpus is what pins behaviour, and it must cover:

| Class | Property |
|---|---|
| Golden positive | valid input produces expected normalised IR |
| Golden negative | invalid input produces the exact code *and* path |
| Permutation invariance | atom array order does not change identity |
| Bond reversal | `{a,b}` ≡ `{b,a}` |
| Round trip | JSON → implementation → JSON preserves semantics |
| Cross-implementation | TypeScript, Python, browser agree |
| Differential | results compared against RDKit |
| Unit metamorphic | equivalent units normalise identically |
| Canonicalisation | equivalent Kekulé forms hash identically |
| Stereo distinction | enantiomers and geometric isomers stay distinct |
| Version migration | older documents migrate predictably |

### 13.2 Chemistry corpus

For the synthesis target, at minimum: water vs hydroperoxyl vs hydrogen
peroxide; benzene in both Kekulé forms; D- vs L-glucose; E- vs Z-2-butene;
formal-charge cases; isotopologues; a supported monoradical; sodium chloride as
multi-component; a hydrate; balanced and deliberately unbalanced equations;
convergent and disconnected routes.

For the CRN target: autocatalysis (`X → X + X`), which must *not* raise a
balance error; a catalytic cycle, which must not raise `route.cycle`; a network
whose deterministic and stochastic semantics visibly diverge at low copy number;
an approximate-majority network; and a bistable switch. Each needs an expected
trajectory, not merely an expected diagnostic — for this target, conformance
includes agreeing on what the network *does*, not only on whether it is
well-formed.

Cross-target cases matter most: the same reaction set read under both profiles
must produce exactly the differences §3.1 predicts, and no others. That test is
what keeps "one shared IR" honest rather than aspirational.

While the v0 profile binds to a named reference implementation, the corpus must
record the reference's version, and any divergence is a specification bug to be
resolved in prose — not a licence to defer to the library indefinitely.

---

## 14. Open questions

1. **Phases and mixtures.** Solutions, suspensions, emulsions have no model.
2. **Polymers.** Repeat units need a bounded representation of an unbounded chain.
3. **Macromolecules.** A protein as an explicit atom list is possible but useless; sequence-level representation is a separate layer, and its absence is why "Bio-Chemical" is not currently an earned description.
4. **Inventory.** Real retrosynthesis needs a modelled catalogue of procurable materials with lots, purity, and cost.
5. **Cross-contamination.** Cleaning between steps is a first-class scheduling constraint with no analogue in the hardware world this IR borrows from.
6. **Provenance.** Literature reference, operator, instrument, and date belong somewhere, and probably not in the chemical layer.
