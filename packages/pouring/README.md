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
networks — targeting both the making of molecules and computing with them.

The IR is language-neutral, with a canonical JSON encoding. This package is the
TypeScript frontend: an embedded DSL that elaborates into that representation,
plus the schema tooling around it. Chemical validation is delegated to an RDKit
reference rather than reimplemented here.



## Status

Early. The chemical IR is specified; the procedure and execution layers are not.
Stereochemistry is unimplemented, so `pouring` cannot yet distinguish D-glucose
from L-glucose, and therefore cannot yet specify a compound unambiguously.



## Syntax

Molecules are graphs, so structure is a netlist — atoms are instances, bonds
connect them by reference:

``` ts
// H₂O — two hydrogens on one oxygen
const water = molecule('water', (m) => {
    const o = m.atom('O');

    m.bond(o, m.atom('H'));
    m.bond(o, m.atom('H'));
});
```

A reaction transforms matter, and is a separate primitive from assembling a
molecule:

``` ts
const acetylation = reaction('acetylation of salicylic acid', {
    reactants: [salicylicAcid, aceticAnhydride],
    products: [aspirin, aceticAcid],
});
```

Quantities carry units, always:

``` ts
const sample = pour(water, 1, 'mol');
```



## What it does not do

`pouring` establishes necessary conditions and never sufficient ones. A
balanced, valence-valid route can still fail in a flask. A program that
validates is not a prediction that a reaction works, and is not authorization
to execute anything physically.



## Documentation

+ [IR specification](https://github.com/daysful/pouring/blob/master/about/notes/ir-specification.md)
+ [Implementation-language decision](https://github.com/daysful/pouring/blob/master/about/notes/implementation-language.md)
