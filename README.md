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
    Bio-Chemical Scripting Language
</h3>



<br />



`pouring` is a high-level language to model and automate the manufacturing of compounds

`pouring` is intended to be

+ used by humans to experiment, test, play, design compounds;
+ used by machines to manufacture compounds;



## Syntax

``` pouring
// the abstract `H_2O`
let waterCompound = react(
    element('H'),
    element('O', 2),
)

// instantiated 1 unit of substance
let waterInstance = pour(
    waterCompound,
    1,
)
```

`react`, `element`, `pour` are primitives


```
// Glucose - https://en.wikipedia.org/wiki/Glucose#/media/File:D-glucose-chain-2D-Fischer.png

let glucoseComponent = react(
    element('H'),
    element('C'),
    react(
        element('O'),
        element('H'),
    ),
)

let glucose = react(
    react(
        element('H'),
        element('C'),
        element('O'),
    ),

    glucoseComponent,

    react(
        react(
            element('H'),
            element('O'),
        ),
        element('C'),
        element('H'),
    ),

    glucoseComponent,

    glucoseComponent,

    react(
        element('C'),
        element('H', 2),
        element('O'),
        element('H'),
    ),
)
```
