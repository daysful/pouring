<p align="center">
    <img src="https://raw.githubusercontent.com/daysful/pouring/master/about/identity/pouring-logo.png" height="250px">
    <br />
    <br />
    <a target="_blank" href="https://github.com/daysful/pouring/blob/master/LICENSE">
        <img src="https://img.shields.io/badge/license-DEL-blue.svg?colorB=1380C3&style=for-the-badge" alt="License: DEL">
    </a>
</p>



<h1 align="center">
    pouring
</h1>


<h3 align="center">
    Chemical Compiler
</h3>



<br />



`pouring` is a typed intermediate representation and compiler toolkit for
chemical structures, reaction equations, material quantities, and reaction
networks.

It targets two things at once, over one shared representation:

+ **making molecules** — synthesis, routes, and eventually procedures compiled
  to laboratory automation;
+ **computing with molecules** — reaction networks as programs, compiled to
  DNA strand displacement.

Both are reaction sets. What differs is the reading laid over them: a synthesis
route is an acyclic plan that terminates, a chemical program is a dynamical
system where cycles are the mechanism. Three rules that look universal hold for
only one of them — mass balance, acyclicity, and whether kinetics matter.

It exists to keep separate the things that are usually run together:

+ what a chemical entity **is**;
+ what physical material is **available**;
+ what transformation is **intended**;
+ what operations are **performed**;
+ how those operations **bind** to a platform.

The IR is language-neutral, with a canonical JSON encoding. TypeScript and
Python eDSLs, structure editors, importers, and a possible future `.pouring`
syntax are all frontends that elaborate into the same checked representation.



## Status

Early. The chemical IR is specified; the procedure and execution layers are not.

There is no hardware, which makes the computational target the nearer one — a
reaction network can be simulated, observed, and checked against what it claims
to compute entirely in software, whereas a synthesis route cannot be verified
without a lab.

**The CRN target runs today.** [`packages/pouring-py`](./packages/pouring-py)
takes a reaction network from definition through validation to both semantics
and back out as a checked result, with no dependencies:

``` bash
cd packages/pouring-py
python3 demo.py
python3 -m unittest discover -s tests
```

The synthesis target covers structure, valence, stereochemistry, units,
balancing, and route reachability. Valence verdicts come from an oracle that
prefers RDKit and falls back to a labelled profile lint — chemistry authority
is delegated by design rather than reimplemented.

Stereochemistry now distinguishes enantiomers: D- and L-glucose share a
formula, a mass, and every bond, and have different identities. What remains
unimplemented is resonance, tautomer, and protonation normalisation, so a
content hash means "the same annotated graph" rather than "the same chemical
entity" — necessary for identity, not sufficient.

See [the IR specification](./about/notes/ir-specification.md) and
[the implementation-language decision](./about/notes/implementation-language.md).



## Syntax

Structures can be written the way chemists write them:

``` python
from pouring import from_smiles

water   = from_smiles('O')
benzene = from_smiles('c1ccccc1')
aspirin = from_smiles('CC(=O)Oc1ccccc1C(=O)O')
```

Aromatic input is resolved to alternating bond orders on the way in, or
refused — the IR has no aromatic bond kind to hide behind.

Underneath, molecules are graphs rather than trees, so structure is a netlist:
atoms are instances and bonds connect them by reference. Because bonds are
declared independently of the order atoms appear, rings close naturally, which
a nested expression tree cannot express:

``` python
from pouring import molecule

def define(m):
    ring = m.atoms('C', 6)
    m.ring(ring, [2, 1, 2, 1, 2, 1])
    m.fill('H')

benzene = molecule('benzene', define)
```

Both paths reach the same molecule — same canonical form, same content hash.

A reaction is a transformation of matter. It is not a way to assemble a
molecule, and the two are different primitives:

``` python
acetylation = Reaction(
    id='acetylation',
    name='acetylation of salicylic acid',
    reactants=[participant('salicylic-acid'), participant('acetic-anhydride')],
    products=[participant('aspirin'), participant('acetic-acid')],
    conditions=Conditions(
        temperature=quantity('90', 'degC'),
        catalyst='sulfuric acid',
    ),
)
```

Quantities carry units, always. One mole and one gram of water are different
objects, and `quantity('90', 'degC')` converts affinely rather than by scaling:

``` python
convert(quantity('90', 'degC'), 'K')     # 363.15 K
convert(quantity('1', 'barg'), 'bar')    # 2.01325 bar absolute
convert(quantity('5', 'g'), 'mL')        # rejected: mass is not volume
```



## What it checks

Valence and connectivity, dimensional consistency, atom and charge balance,
route reachability from declared starting materials.



## What it does not

`pouring` establishes **necessary** conditions and never sufficient ones. A
balanced, valence-valid, well-typed route can still fail in a flask — chemistry
does not compose the way logic gates do, and protecting groups exist precisely
because functional groups interfere non-locally.

A program that validates is not a prediction that a reaction works, and is not
authorization to execute anything physically.
