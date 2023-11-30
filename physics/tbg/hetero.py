#! /usr/bin/env python
"""
Example for constructing twisted bilayer graphene.

References:
[1] https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.256802
[2] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.125413
"""

from math import acos

import numpy as np
from ase.io import read, write

from utils import frac2cart, twist_cell, make_layer


def calc_twist_angle(i: int) -> float:
    """
    Calculate twisting angle according to ref. [1].

    :param i: twist index
    :return: twisting angle in RADIANs, NOT degrees
    """
    cos_ang = (3 * i ** 2 + 3 * i + 0.5) / (3 * i ** 2 + 3 * i + 1)
    return acos(cos_ang)


def calc_twist_angle2(n: int, m: int) -> float:
    """
    Calculate twisting angle according to ref. [2].

    :param n: twist index n
    :param m: twist index m
    :return: twisting angle in RADIANs, NOT degrees
    """
    cos_ang = (n**2 + 4 * n * m + m**2) / (2 * (n**2 + n * m + m**2))
    return acos(cos_ang)


def calc_hetero_lattice(i: int) -> np.ndarray:
    """
    Calculate fractional coordinates of lattice vectors of hetero-structure
    according to ref. [1].

    :param i: parameter controlling the twisting angle
    :return: (3, 3) float64 array
    """
    hetero_lattice = np.array([[i, i + 1, 0],
                               [-(i + 1), 2 * i + 1, 0],
                               [0, 0, 1]])
    return hetero_lattice


def calc_hetero_lattice2(n: int, m: int) -> np.ndarray:
    """
    Calculate fractional coordinates of lattice vectors of hetero-structure
    according to ref. [2].

    :param n: parameter controlling the twisting angle
    :param m: parameter controlling the twisting angle
    :return: (3, 3) float64 array
    """
    hetero_lattice = np.array([[n, m, 0],
                               [-m, n + m, 0],
                               [0, 0, 1]])
    return hetero_lattice


def main():
    # Evaluate twisting angle and hetero lattice vectors
    i = 5
    angle = calc_twist_angle(i)
    hetero_lattice = calc_hetero_lattice(i)

    # Load primitive cells of each layer
    cell_fixed = read("POSCAR", index=-1)
    cell_twisted = cell_fixed.copy()

    # Twist primitive cell of top layer
    twist_cell(cell_twisted, angle=angle, shift=3.349)

    # Evaluate Cartesian coordinates of lattice vectors of hetero-structure
    hetero_lattice = frac2cart(cell_fixed.cell, hetero_lattice)
    layer_fixed = make_layer(cell_fixed, hetero_lattice)
    layer_twisted = make_layer(cell_twisted, hetero_lattice)

    # Merge layers
    tbg = layer_fixed + layer_twisted

    # Visualization
    tbg *= (3, 3, 1)
    vasp_args = {"direct": True, "sort": True}
    write("POSCAR.tbg", tbg, **vasp_args)


if __name__ == "__main__":
    main()
