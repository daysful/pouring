// #region module
/**
 * The `pouring` intermediate representation.
 *
 * Three layers, deliberately separated:
 *
 *  1. structure  — what a compound IS       (Fragment: atoms + bonds + ports)
 *  2. process    — how a compound is MADE   (Reaction, Route)
 *  3. quantity   — how much of it EXISTS    (Quantity, Pour)
 *
 * Everything here is plain serialisable data. Builders produce it,
 * validators check it, backends consume it. No behaviour lives on these types.
 */



// #region structure

export type ElementSymbol = string;

export type AtomId = string;
export type BondId = string;
export type PortId = string;


/**
 * Bond order in whole units. Aromatic systems are stored in Kekule form
 * (alternating 1/2) so that valence arithmetic stays exact; aromaticity is a
 * perception pass over the graph, not a bond kind.
 */
export type BondOrder = 1 | 2 | 3;


export interface Atom {
    id: AtomId;
    element: ElementSymbol;
    /** Formal charge. Shifts the effective valence. */
    charge: number;
    /** Mass number, when a specific isotope is intended. */
    isotope?: number;
    /** Human-facing name, used in diagnostics only. */
    label?: string;
}


export interface Bond {
    id: BondId;
    from: AtomId;
    to: AtomId;
    order: BondOrder;
}


/**
 * An unsatisfied valence deliberately left open, so a fragment can be wired
 * into a larger structure later. The HDL parallel is exact: a port is an
 * unconnected pin, `width` is its bus width, and elaboration connects it.
 *
 * An open valence that is NOT declared as a port is an error (a floating net).
 * An open valence that IS declared is either an attachment point or, on a
 * complete molecule, a radical.
 */
export interface Port {
    id: PortId;
    name: string;
    atom: AtomId;
    width: BondOrder;
}


/**
 * A molecular graph. A Fragment with no open ports is a complete molecule —
 * the elaborated, manufacturable thing. A Fragment with ports is a reusable
 * module: a functional group, a protecting group, a polymer repeat unit.
 */
export interface Fragment {
    kind: 'fragment';
    name: string;
    atoms: Atom[];
    bonds: Bond[];
    ports: Port[];
}

// #endregion structure



// #region quantity

/**
 * Units carry their dimension so the checker can reject `5 grams + 2 kelvin`
 * before any of it reaches a machine.
 */
export type Dimension =
    | 'amount'
    | 'mass'
    | 'volume'
    | 'temperature'
    | 'time'
    | 'concentration';

export type Unit =
    | 'mol' | 'mmol'
    | 'g' | 'mg' | 'kg'
    | 'L' | 'mL'
    | 'K' | 'C'
    | 's' | 'min' | 'h'
    | 'M';

export interface Quantity {
    value: number;
    unit: Unit;
}


/** A concrete physical instance of a compound — the result of `pour`. */
export interface Pour {
    kind: 'pour';
    species: string;
    amount: Quantity;
}

// #endregion quantity



// #region process

export interface Stoichiometry {
    /** Name of a Fragment with no open ports. */
    species: string;
    coefficient: number;
}


export interface Conditions {
    temperature?: Quantity;
    time?: Quantity;
    pressure?: Quantity;
    solvent?: string;
    catalyst?: string;
}


/**
 * A transformation of matter: the process primitive.
 *
 * Note this is what `react` means here — reactants become products under
 * conditions. Composing atoms into a molecule is NOT a reaction; that is
 * structure, and it lives in Fragment.
 */
export interface Reaction {
    kind: 'reaction';
    name: string;
    reactants: Stoichiometry[];
    products: Stoichiometry[];
    conditions?: Conditions;
    /** Fractional yield in [0, 1]. Absent means unknown, not quantitative. */
    yield?: number;
}


/**
 * A synthesis plan: species in scope plus the reactions that consume and
 * produce them. The reaction list forms a DAG once resolved — this is the
 * netlist of the process layer.
 */
export interface Route {
    kind: 'route';
    name: string;
    target: string;
    species: Fragment[];
    reactions: Reaction[];
}

// #endregion process



// #region diagnostics

export type Severity = 'error' | 'warning';

export interface Diagnostic {
    severity: Severity;
    /** Stable machine-readable code, e.g. 'valence.floating'. */
    code: string;
    message: string;
    /** Whatever the diagnostic is anchored to: atom id, reaction name, ... */
    subject?: string;
}

// #endregion diagnostics
// #endregion module
