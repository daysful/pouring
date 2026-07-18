// #region module
import {
    ElementSymbol,
} from './types';



// #region data

export interface ElementData {
    symbol: ElementSymbol;
    number: number;
    /** Standard atomic weight, IUPAC 2021 abridged. */
    mass: number;
    /** Common valences, in the order they should be preferred when filling. */
    valences: number[];
}


/**
 * Deliberately small. `pouring` should cover organic and simple inorganic
 * chemistry correctly before it covers everything badly.
 */
export const ELEMENTS: Record<ElementSymbol, ElementData> = {
    H:  { symbol: 'H',  number: 1,  mass: 1.008,   valences: [1] },
    He: { symbol: 'He', number: 2,  mass: 4.0026,  valences: [0] },
    B:  { symbol: 'B',  number: 5,  mass: 10.81,   valences: [3] },
    C:  { symbol: 'C',  number: 6,  mass: 12.011,  valences: [4] },
    N:  { symbol: 'N',  number: 7,  mass: 14.007,  valences: [3, 5] },
    O:  { symbol: 'O',  number: 8,  mass: 15.999,  valences: [2] },
    F:  { symbol: 'F',  number: 9,  mass: 18.998,  valences: [1] },
    Ne: { symbol: 'Ne', number: 10, mass: 20.180,  valences: [0] },
    Na: { symbol: 'Na', number: 11, mass: 22.990,  valences: [1] },
    Mg: { symbol: 'Mg', number: 12, mass: 24.305,  valences: [2] },
    Al: { symbol: 'Al', number: 13, mass: 26.982,  valences: [3] },
    Si: { symbol: 'Si', number: 14, mass: 28.085,  valences: [4] },
    P:  { symbol: 'P',  number: 15, mass: 30.974,  valences: [3, 5] },
    S:  { symbol: 'S',  number: 16, mass: 32.06,   valences: [2, 4, 6] },
    Cl: { symbol: 'Cl', number: 17, mass: 35.45,   valences: [1] },
    Ar: { symbol: 'Ar', number: 18, mass: 39.948,  valences: [0] },
    K:  { symbol: 'K',  number: 19, mass: 39.098,  valences: [1] },
    Ca: { symbol: 'Ca', number: 20, mass: 40.078,  valences: [2] },
    Fe: { symbol: 'Fe', number: 26, mass: 55.845,  valences: [2, 3] },
    Cu: { symbol: 'Cu', number: 29, mass: 63.546,  valences: [1, 2] },
    Zn: { symbol: 'Zn', number: 30, mass: 65.38,   valences: [2] },
    Br: { symbol: 'Br', number: 35, mass: 79.904,  valences: [1] },
    I:  { symbol: 'I',  number: 53, mass: 126.90,  valences: [1] },
};

// #endregion data



// #region access

export const isKnownElement = (
    symbol: ElementSymbol,
): boolean => Object.prototype.hasOwnProperty.call(ELEMENTS, symbol);


export const elementData = (
    symbol: ElementSymbol,
): ElementData => {
    const data = ELEMENTS[symbol];

    if (!data) {
        throw new Error(`[pouring] unknown element '${symbol}'`);
    }

    return data;
};


/**
 * Valences available to an atom, adjusted for formal charge.
 *
 * A charge changes how many bonds an atom can carry: ammonium N+ takes 4
 * rather than 3, hydroxide O- takes 1 rather than 2. The sign convention
 * differs across the periodic table; this covers the main-group organic cases
 * `pouring` targets today.
 */
export const availableValences = (
    symbol: ElementSymbol,
    charge: number,
): number[] => {
    const { valences } = elementData(symbol);

    if (charge === 0) {
        return valences;
    }

    return valences
        .map((valence) => valence + charge)
        .filter((valence) => valence >= 0);
};

// #endregion access
// #endregion module
