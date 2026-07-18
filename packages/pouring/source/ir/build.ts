// #region module
import {
    Atom,
    AtomId,
    Bond,
    BondOrder,
    Conditions,
    ElementSymbol,
    Fragment,
    Port,
    Pour,
    Quantity,
    Reaction,
    Route,
    Stoichiometry,
    Unit,
} from './types';

import {
    availableValences,
    elementData,
} from './elements';



// #region structure

export interface AtomOptions {
    charge?: number;
    isotope?: number;
    label?: string;
}


/**
 * Construction context for a single compound.
 *
 * Atoms are instances, bonds are nets, ports are unconnected pins. Because
 * nets are named independently of the order atoms are declared in, cyclic
 * structures fall out for free — a ring is just a net that closes back.
 */
export class Compound {
    private atomsList: Atom[] = [];
    private bondsList: Bond[] = [];
    private portsList: Port[] = [];
    private counter = 0;

    constructor(
        public readonly name: string,
    ) {}


    private nextId(
        prefix: string,
    ): string {
        return `${prefix}${this.counter++}`;
    }


    atom(
        element: ElementSymbol,
        options: AtomOptions = {},
    ): AtomId {
        elementData(element);

        const atom: Atom = {
            id: this.nextId('a'),
            element,
            charge: options.charge ?? 0,
        };

        if (options.isotope !== undefined) {
            atom.isotope = options.isotope;
        }

        if (options.label !== undefined) {
            atom.label = options.label;
        }

        this.atomsList.push(atom);

        return atom.id;
    }


    atoms(
        element: ElementSymbol,
        count: number,
        options: AtomOptions = {},
    ): AtomId[] {
        return Array.from(
            { length: count },
            () => this.atom(element, options),
        );
    }


    bond(
        from: AtomId,
        to: AtomId,
        order: BondOrder = 1,
    ): Bond {
        if (from === to) {
            throw new Error(
                `[pouring] '${this.name}': cannot bond atom '${from}' to itself`,
            );
        }

        this.assertAtom(from);
        this.assertAtom(to);

        const bond: Bond = {
            id: this.nextId('b'),
            from,
            to,
            order,
        };

        this.bondsList.push(bond);

        return bond;
    }


    /** Bonds atoms in an open sequence: n atoms, n-1 bonds. */
    chain(
        atoms: AtomId[],
        orders?: BondOrder[],
    ): Bond[] {
        return atoms.slice(0, -1).map(
            (atom, index) => this.bond(
                atom,
                atoms[index + 1],
                orders?.[index] ?? 1,
            ),
        );
    }


    /**
     * Bonds atoms in a closed cycle: n atoms, n bonds.
     *
     * This is the operation a nested-call tree cannot express, and the reason
     * the IR is a netlist rather than an expression tree.
     */
    ring(
        atoms: AtomId[],
        orders?: BondOrder[],
    ): Bond[] {
        if (atoms.length < 3) {
            throw new Error(
                `[pouring] '${this.name}': a ring needs at least 3 atoms`,
            );
        }

        return atoms.map(
            (atom, index) => this.bond(
                atom,
                atoms[(index + 1) % atoms.length],
                orders?.[index] ?? 1,
            ),
        );
    }


    /** Declares an intentionally open valence: an attachment point. */
    port(
        atom: AtomId,
        name: string,
        width: BondOrder = 1,
    ): Port {
        this.assertAtom(atom);

        if (this.portsList.some((existing) => existing.name === name)) {
            throw new Error(
                `[pouring] '${this.name}': duplicate port '${name}'`,
            );
        }

        const port: Port = {
            id: this.nextId('p'),
            name,
            atom,
            width,
        };

        this.portsList.push(port);

        return port;
    }


    /**
     * Elaborates a fragment into this compound, wiring its ports to local
     * atoms — module instantiation, in the HDL sense. Interior atoms are
     * copied with fresh ids so the same fragment can be used many times.
     */
    instance(
        fragment: Fragment,
        connections: Record<string, AtomId>,
    ): AtomId[] {
        const remap = new Map<AtomId, AtomId>();

        for (const atom of fragment.atoms) {
            const copy: Atom = {
                ...atom,
                id: this.nextId('a'),
            };

            this.atomsList.push(copy);
            remap.set(atom.id, copy.id);
        }

        for (const bond of fragment.bonds) {
            this.bondsList.push({
                ...bond,
                id: this.nextId('b'),
                from: remap.get(bond.from)!,
                to: remap.get(bond.to)!,
            });
        }

        for (const port of fragment.ports) {
            const target = connections[port.name];

            if (target === undefined) {
                throw new Error(
                    `[pouring] '${this.name}': port '${port.name}' of `
                    + `'${fragment.name}' left unconnected`,
                );
            }

            this.bond(remap.get(port.atom)!, target, port.width);
        }

        for (const name of Object.keys(connections)) {
            if (!fragment.ports.some((port) => port.name === name)) {
                throw new Error(
                    `[pouring] '${this.name}': '${fragment.name}' has no `
                    + `port '${name}'`,
                );
            }
        }

        return Array.from(remap.values());
    }


    /**
     * Saturates every unfilled valence with the given element, hydrogen by
     * default. This is the implicit-hydrogen convention every chemist writes
     * in — the IR stays fully explicit, the source stays readable.
     *
     * Atoms whose valence is already met, exceeded, or spoken for by a port
     * are left alone; over-filled atoms are reported by the validator rather
     * than silently repaired here.
     */
    fill(
        element: ElementSymbol = 'H',
    ): void {
        const existing = [...this.atomsList];

        for (const atom of existing) {
            const used = this.usedValence(atom.id);
            const valences = availableValences(atom.element, atom.charge);
            const target = valences.find((valence) => valence >= used);

            if (target === undefined) {
                continue;
            }

            for (let index = 0; index < target - used; index++) {
                this.bond(atom.id, this.atom(element));
            }
        }
    }


    private usedValence(
        atom: AtomId,
    ): number {
        const bonded = this.bondsList
            .filter((bond) => bond.from === atom || bond.to === atom)
            .reduce((total, bond) => total + bond.order, 0);

        const ported = this.portsList
            .filter((port) => port.atom === atom)
            .reduce((total, port) => total + port.width, 0);

        return bonded + ported;
    }


    private assertAtom(
        atom: AtomId,
    ): void {
        if (!this.atomsList.some((existing) => existing.id === atom)) {
            throw new Error(
                `[pouring] '${this.name}': unknown atom '${atom}'`,
            );
        }
    }


    freeze(): Fragment {
        return {
            kind: 'fragment',
            name: this.name,
            atoms: [...this.atomsList],
            bonds: [...this.bondsList],
            ports: [...this.portsList],
        };
    }
}


/**
 * Declares a compound.
 *
 * A compound with no open ports is a complete molecule. A compound with ports
 * is a reusable module — a functional group, a protecting group, a repeat unit.
 */
export const compound = (
    name: string,
    define: (compound: Compound) => void,
): Fragment => {
    const building = new Compound(name);

    define(building);

    return building.freeze();
};

// #endregion structure



// #region process

export type Reagent = Fragment | [number, Fragment];


const stoichiometry = (
    reagent: Reagent,
): Stoichiometry => Array.isArray(reagent)
    ? { species: reagent[1].name, coefficient: reagent[0] }
    : { species: reagent.name, coefficient: 1 };


export interface ReactionSpecification {
    reactants: Reagent[];
    products: Reagent[];
    conditions?: Conditions;
    yield?: number;
}


/**
 * Declares a reaction: reactants become products under conditions.
 *
 * `react` is the process primitive. Assembling atoms into a molecule is not a
 * reaction — that is structure, and it belongs to `compound`.
 */
export const react = (
    name: string,
    specification: ReactionSpecification,
): Reaction => {
    const reaction: Reaction = {
        kind: 'reaction',
        name,
        reactants: specification.reactants.map(stoichiometry),
        products: specification.products.map(stoichiometry),
    };

    if (specification.conditions !== undefined) {
        reaction.conditions = specification.conditions;
    }

    if (specification.yield !== undefined) {
        reaction.yield = specification.yield;
    }

    return reaction;
};


export const route = (
    name: string,
    target: Fragment,
    species: Fragment[],
    reactions: Reaction[],
): Route => ({
    kind: 'route',
    name,
    target: target.name,
    species,
    reactions,
});

// #endregion process



// #region quantity

/**
 * Instantiates a physical amount of a compound.
 *
 * The unit is required. `pour(water, 1)` is not a quantity of anything —
 * one mole and one gram of water are different objects, and a language that
 * drives real hardware cannot be vague about which.
 */
export const pour = (
    species: Fragment,
    value: number,
    unit: Unit,
): Pour => ({
    kind: 'pour',
    species: species.name,
    amount: { value, unit },
});


export const quantity = (
    value: number,
    unit: Unit,
): Quantity => ({ value, unit });

// #endregion quantity
// #endregion module
