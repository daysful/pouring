# pouring

## biochemical scripting language


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
