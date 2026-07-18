# Implementation language

**Status:** proposed, revised
**Decision:** TypeScript for frontend and IR tooling. Python plus RDKit as the
reference chemistry implementation, starting now. Rust deferred behind
evidence-based triggers.

> **Revision note.** The previous draft postponed Python to "an eventual
> authoring layer" and treated the TypeScript implementation as the chemistry
> authority. That was wrong, and §3 explains why. The Rust triggers have also
> been rewritten from feature-based to evidence-based.

---

## 1. The question is smaller than it looks

Because [the IR](./ir-specification.md) is specified as language-neutral data
with a JSON encoding, no implementation owns it. Whatever language writes the
first validator is writing *an* implementation, not *the* implementation, and a
second one in another language stays interoperable by construction.

That reframes the decision from "what is `pouring` written in" — near
irreversible — to "what do we prototype in", which is cheap to change. Making
the spec language-neutral first is what buys that option.

## 2. What the system needs

1. **Graph work** — molecular graphs, ring perception, canonicalisation
2. **Substructure matching** — subgraph isomorphism; NP-hard
3. **Search** — retrosynthesis over a large transform space
4. **Exact arithmetic** — rational stoichiometry, nullspace balancing
5. **Chemistry semantics** — valence, aromaticity perception, sanitisation, stereo
6. **Compiler ergonomics** — sum types, exhaustive matching, immutable transforms
7. **Backend emission** — XDL (XML), Autoprotocol (JSON), ORD interchange
8. **Numerical simulation** — mass-action ODE integration and stochastic
   (Gillespie) simulation of reaction networks

Items 1, 4, 6, 7 are comfortable anywhere. Items 2 and 3 are where performance
bites, and neither is near-term. **Item 5 is the one that changes the decision**,
and the previous draft missed it by filing everything under performance.

Item 8 arrived with the CRN target and points the same way. Numerical ODE
integration and stochastic simulation are exactly what SciPy exists for, and
have no comparable JavaScript story. More importantly, with no hardware
available, **simulation is the only verification path that actually closes** —
a reaction network can be run and checked against what it claims to compute,
while a synthesis route cannot be verified without a lab. That makes the
numerics load-bearing rather than incidental, and it lands in Python's column.

## 3. The correction: chemistry is a correctness problem, not a performance one

The previous draft reasoned entirely about speed, concluded nothing near-term
was compute-bound, and therefore chose TypeScript for everything. The reasoning
about speed was right. The conclusion was wrong, because the hard part of
items 1 and 5 was never speed.

Aromaticity perception, sanitisation, valence special cases, stereo perception,
canonicalisation, and SMARTS semantics are **decades of accumulated edge cases**.
RDKit's sanitisation alone runs several context-dependent passes — valence
checking, organometallic cleanup, Kekulisation, aromaticity perception, radical
assignment, chirality cleanup — and its documented valence handling is
substantially more nuanced than any table one would write from memory.

The evidence is already in this repository. The probe's `source/ir/elements.ts`
contains a hand-written valence table asserting `S: [2,4,6]` and `N: [3,5]`.
That is a crude approximation of something with real subtlety, written without
a reference, and it would have silently produced wrong answers rather than
failing loudly. The specification now concedes the same point directly: its
valence rule is "not specified precisely enough for two independent
implementations to agree" and binds to a named reference.

**Reimplementing chemistry in TypeScript is therefore not a shortcut. It is the
most expensive possible way to be subtly wrong**, and it is expensive in the
currency that matters least to notice and most to fix.

## 4. Candidates

**TypeScript.** The scaffold exists. Discriminated unions model IR nodes well,
and the eDSL pattern — a library that builds a netlist when executed — is proven
by Chisel and Amaranth. Compiles to the browser, which makes an in-browser
playground possible; for a language nobody has heard of, "try it without
installing anything" is a real adoption lever. Weak numerics, slow for search,
and — decisively — **no credible path to owning chemistry semantics**.

The usual objection, "no chemistry ecosystem", is subtler than it sounds. RDKit
ships an official WASM build, so JavaScript can *call* RDKit. That is an argument
for TypeScript as a frontend that delegates, not as the place where chemistry is
defined.

**Python.** The native language of computational chemistry: RDKit, OpenBabel,
PySCF, ASKCOS, AiZynthFinder. Amaranth proves the eDSL pattern here too, and
Jupyter is the medium chemists work in. Direct, complete, testable access to the
reference toolkit. Slow, and packaging is a persistent tax.

**Rust.** The strongest endgame. Sum types and exhaustive matching suit compiler
work; speed suits search and substructure matching. The decisive property is
distribution: one core to WASM, to Python via PyO3, to Node via napi-rs — the
architecture polars, ruff, and tokenizers converged on. Costs the most to write,
and iterating an unsettled design in Rust is slow.

**Haskell / OCaml.** Best in class for IR transforms, and Clash shows an HDL can
live here. Ruled out on collaboration and ecosystem grounds, not technical ones.

## 5. Decision

**Split by what each language is authoritative for.**

| TypeScript | Python + RDKit |
|---|---|
| eDSL and elaboration | reference chemical validation |
| JSON Schema tooling | normalisation and sanitisation |
| diagnostics presentation | identifiers — InChI, canonical SMILES |
| editor, language server | formula and molar mass |
| browser playground | stereochemistry perception |
| chemistry-*independent* IR transforms | canonicalisation experiments |
| backend adapters | ORD and dataset interoperability |
| trajectory visualisation | CRN simulation — mass-action ODE, Gillespie |

The dividing line is not performance and not seniority. It is **which side owns
a claim about chemistry**. TypeScript may call RDKit's WASM build in the browser;
conformance tests then compare those results against the Python reference, and
divergence is a bug in the specification rather than a tolerance to live with.

**Move the core to Rust on evidence, not on features.** The previous triggers
("retrosynthesis exists") were feature-based and would have fired prematurely —
retrosynthesis can perfectly well remain a Python service or delegate to an
existing engine. Replace them with:

| Trigger | Evidence required |
|---|---|
| Performance | representative profiling shows the core dominates runtime |
| Distribution | users need one local binary across platforms |
| Divergence | TS, Python, and browser implementations repeatedly disagree |
| Runtime reliability | a long-running execution engine needs strict resource and concurrency control |
| Stable semantics | IR v1 and conformance behaviour are settled |
| Maintenance economics | one Rust core is demonstrably cheaper than several implementations |

Until one is *measured*, Rust buys performance nobody is waiting on and costs
iteration speed the design needs.

## 6. Consequences

**The implementations are replaceable and conformance-bound.** The previous
draft called the TypeScript work "disposable", which invites bad engineering and
undersells it — the eDSL, schema tooling, editor, tests, and examples stay
valuable even if a core is later written in Rust. Held to the conformance suite,
any implementation may be replaced; none may be sloppy.

**The spec stays ahead of every implementation.** Where they disagree, the spec
wins and the code is the bug.

**The conformance suite must exist before the second implementation does.** It is
what makes "many frontends, many backends" true rather than aspirational, and
with two languages now in scope from the start, it is needed immediately rather
than eventually. Requirements are in the specification, §13.

**The existing probe's chemistry is contra-indicated.** `source/ir/elements.ts`
should not grow; its valence table is exactly the kind of assertion that belongs
behind an RDKit reference. The probe remains useful for IR shape and eDSL
ergonomics, which is what it was for.

## 7. Open

**Chemists are not one audience.** "Chemists use Python" is directionally right
for computational chemistry and wrong for bench synthetic chemists, who are
better served by graphical editors, ELN import, or assisted authoring. Python is
necessary, not sufficient, and it is not the only user-facing path.

**Whether `.pouring` ever gets a real parser** is a separate and later decision.
Chisel and Amaranth never grew custom syntax and are not obviously worse for it.
A standalone syntax earns its cost when non-programmers must author, when static
analysis must happen without executing arbitrary code, or when an archival
textual form is required — an eDSL source file can read the network, use the
clock, and depend on unpinned packages, so it is not by itself a reproducible
scientific record. The canonical artefact is the normalised IR plus schema
version, chemistry profile, compiler version, dependency lock, and source hash.
