#! /usr/bin/env python
"""
Example for constructing twisted bilayer graphene.

References:
[1] https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.256802
[2] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.125413
"""

import math

import numpy as np
from numpy.linalg import norm

import tbplas as tb


def calc_twist_angle(i: int) -> float:
    """
    Calculate twisting angle according to ref. [1].

    :param i: parameter controlling the twisting angle
    :return: twisting angle in RADIANs, NOT degrees
    """
    cos_ang = (3 * i ** 2 + 3 * i + 0.5) / (3 * i ** 2 + 3 * i + 1)
    return math.acos(cos_ang)


def calc_twist_angle2(n: int, m: int) -> float:
    """
    Calculate twisting angle according to ref. [2].

    :param n: parameter controlling the twisting angle
    :param m: parameter controlling the twisting angle
    :return: twisting angle in RADIANs, NOT degrees
    """
    cos_ang = (n**2 + 4 * n * m + m**2) / (2 * (n**2 + n * m + m**2))
    return math.acos(cos_ang)


def calc_hetero_lattice(i: int, prim_cell_fixed: tb.PrimitiveCell) -> np.ndarray:
    """
    Calculate Cartesian coordinates of lattice vectors of hetero-structure
    according to ref. [1].

    :param i: parameter controlling the twisting angle
    :param prim_cell_fixed: primitive cell of fixed layer
    :return: (3, 3) float64 array, Cartesian coordinates of hetero-structure
        lattice vectors in NANOMETER
    """
    hetero_lattice = np.array([[i, i + 1, 0],
                               [-(i + 1), 2 * i + 1, 0],
                               [0, 0, 1]])
    hetero_lattice = tb.frac2cart(prim_cell_fixed.lat_vec, hetero_lattice)
    return hetero_lattice


def calc_hetero_lattice2(n: int, m: int,
                         prim_cell_fixed: tb.PrimitiveCell) -> np.ndarray:
    """
    Calculate Cartesian coordinates of lattice vectors of hetero-structure
    according to ref. [2].

    :param n: parameter controlling the twisting angle
    :param m: parameter controlling the twisting angle
    :param prim_cell_fixed: primitive cell of fixed layer
    :return: (3, 3) float64 array, Cartesian coordinates of hetero-structure
        lattice vectors in NANOMETER
    """
    hetero_lattice = np.array([[n, m, 0],
                               [-m, n + m, 0],
                               [0, 0, 1]])
    hetero_lattice = tb.frac2cart(prim_cell_fixed.lat_vec, hetero_lattice)
    return hetero_lattice


def calc_hop(rij: np.ndarray) -> float:
    """
    Calculate hopping parameter according to Slater-Koster relation.
    See ref. [2] for the formulae.

    :param rij: (3,) array, displacement vector between two orbitals in NM
    :return: hopping parameter in eV
    """
    a0 = 0.1418
    a1 = 0.3349
    r_c = 0.6140
    l_c = 0.0265
    gamma0 = 2.7
    gamma1 = 0.48
    decay = 22.18
    q_pi = decay * a0
    q_sigma = decay * a1
    dr = norm(rij).item()
    n = rij.item(2) / dr
    v_pp_pi = - gamma0 * math.exp(q_pi * (1 - dr / a0))
    v_pp_sigma = gamma1 * math.exp(q_sigma * (1 - dr / a1))
    fc = 1 / (1 + math.exp((dr - r_c) / l_c))
    hop = (n**2 * v_pp_sigma + (1 - n**2) * v_pp_pi) * fc
    return hop


def extend_hop(prim_cell: tb.PrimitiveCell, max_distance: float = 0.75) -> None:
    """
    Extend the hopping terms in primitive cell up to cutoff distance.

    :param prim_cell: primitive cell to extend
    :param max_distance: cutoff distance in NM
    :return: None. Incoming primitive cell is modified
    """
    neighbors = tb.find_neighbors(prim_cell, a_max=1, b_max=1,
                                  max_distance=max_distance)
    for term in neighbors:
        i, j = term.pair
        prim_cell.add_hopping(term.rn, i, j, calc_hop(term.rij))


def make_tbg(i: int = 1, shift: float = 0.3349) -> tb.PrimitiveCell:
    """
    Make primitive cell of TBG.

    :param i: integer controlling the twisting angle
    :param shift: interlayer distance in nm
    :return: primitive cell of TBG.
    """
    # Evaluate twisting angle.
    angle = calc_twist_angle(i)

    # Build layers
    prim_cell_fixed = tb.make_graphene_diamond()
    prim_cell_twisted = tb.make_graphene_diamond()

    # Shift and twist the primitive cell of top layer
    tb.spiral_prim_cell(prim_cell_twisted, angle=angle, shift=shift)

    # Shift the on-site energies of both primitive cells by -0.78 eV according
    # to the appendix of Ref 2.
    for cell in (prim_cell_fixed, prim_cell_twisted):
        for orb_i in range(cell.num_orb):
            cell.set_orbital(orb_i, energy=-0.78)

    # Evaluate coordinates of lattice vectors of hetero-structure.
    hetero_lattice = calc_hetero_lattice(i, prim_cell_fixed)

    # Reshape primtive cells to produce layers
    layer_fixed = tb.make_hetero_layer(prim_cell_fixed, hetero_lattice)
    layer_twisted = tb.make_hetero_layer(prim_cell_twisted, hetero_lattice)

    # Merge layers and extend hopping terms
    merged_cell = tb.merge_prim_cell(layer_fixed, layer_twisted)
    extend_hop(merged_cell, max_distance=0.5)
    angle = -math.atan(hetero_lattice[0, 1] / hetero_lattice[0, 0])
    tb.spiral_prim_cell(merged_cell, angle=angle)
    return merged_cell


def main(i: int = 25) -> None:
    tbg = make_tbg(i)
    super_cell = tb.SuperCell(tbg, dim=(3, 3, 1), pbc=(True, True, False))
    sample = tb.Sample(super_cell)
    sample.plot(with_conj=False, with_orbitals=False, hop_as_arrows=False,
                hop_eng_cutoff=0.5)


if __name__ == "__main__":
    main()
