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

Early, and the honest summary is short: the chemical IR is specified, the
procedure and execution layers are not, and there is no working implementation
yet.

Stereochemistry is unimplemented. Until it lands, `pouring` cannot distinguish
D-glucose from L-glucose — they have identical atoms, bonds, and formula — and
therefore cannot yet specify a compound unambiguously.

There is no hardware. That makes the computational target the nearer one: a
reaction network can be simulated, observed, and checked against what it claims
to compute entirely in software, whereas a synthesis route cannot be verified
without a lab.

See [the IR specification](./about/notes/ir-specification.md) and
[the implementation-language decision](./about/notes/implementation-language.md).



## Syntax

Molecules are graphs, not trees, so structure is written as a netlist: atoms are
instances, and bonds connect them by reference.

``` ts
// H₂O — two hydrogens on one oxygen
const water = molecule('water', (m) => {
    const o = m.atom('O');

    m.bond(o, m.atom('H'));
    m.bond(o, m.atom('H'));
});
```

Because bonds are declared independently of the order atoms appear, rings close
naturally — which a nested expression tree cannot express:

``` ts
const benzene = molecule('benzene', (m) => {
    const ring = m.atoms('C', 6);

    m.ring(ring, [2, 1, 2, 1, 2, 1]);
    m.fill('H');
});
```

A reaction is a transformation of matter. It is not a way to assemble a
molecule, and the two are different primitives:

``` ts
const acetylation = reaction('acetylation of salicylic acid', {
    reactants: [salicylicAcid, aceticAnhydride],
    products: [aspirin, aceticAcid],
    conditions: {
        temperature: quantity(90, 'degC'),
        catalyst: sulfuricAcid,
    },
});
```

Quantities carry units, always. One mole and one gram of water are different
objects:

``` ts
const sample = pour(water, 1, 'mol');
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
