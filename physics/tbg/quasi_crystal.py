#! /usr/bin/env python
"""
Example for constructing quasi-crystal at primitive cell and sample levels.
"""

import math
from typing import Tuple

import numpy as np
from numpy.linalg import norm
from ase import Atoms
from ase.io import read, write

from utils import rotate_coord, twist_cell


def cutoff(prim_cell: Atoms, center: np.ndarray, radius: float = 30.0) -> None:
    """
    Cutoff primitive cell up to given radius with respect to given center.

    :param prim_cell: primitive cell to cut
    :param center: Cartesian coordinate of center in angstrom
    :param radius: cutoff radius in angstrom
    :return: None. Incoming primitive cell is modified.
    """
    idx_remove = []
    orb_pos = prim_cell.get_positions()
    for i, pos in enumerate(orb_pos):
        if norm(pos[:2] - center[:2]) > radius:
            idx_remove.append(i)
    del prim_cell[idx_remove]


def make_quasi_crystal(prim_cell: Atoms,
                       dim: Tuple[int, int, int],
                       angle: float,
                       center: np.ndarray,
                       radius: float = 30.0,
                       shift: float = 3.0,
                       algo: int = 1):
    """
    Create quasi-crystal.

    :param prim_cell: primitive cell from which the quasi-crystal is built
    :param dim: dimension of the extended primitive cell
    :param angle: twisting angle in RADIAN
    :param center: fractional coordinate of the twisting center in the primitive
        cell
    :param radius: radius of quasi-crystal in angstrom
    :param shift: distance of shift along z-axis in angstrom
    :param algo: algorithm to build the quasi-crystal
    :return: None
    """
    # Get the Cartesian coordinate of rotation center
    center = np.array([dim[0]//2, dim[1]//2, 0]) + center
    center = np.matmul(center, prim_cell.cell)

    # Build fixed and twisted layers
    layer_fixed = prim_cell * dim
    layer_twisted = prim_cell * dim

    if algo == 1:
        # Rotate and shift twisted layer
        twist_cell(layer_twisted, angle=angle, center=center, shift=shift)

        # Remove unnecessary orbitals
        cutoff(layer_fixed, center=center, radius=radius)
        cutoff(layer_twisted, center=center, radius=radius)

        # Reset the lattice of twisted layer
        layer_twisted.set_cell(layer_fixed.cell, scale_atoms=False)
        write("POSCAR.qs", layer_fixed + layer_twisted)
    else:
        # Remove unnecessary orbitals
        cutoff(layer_fixed, center=center, radius=radius)
        cutoff(layer_twisted, center=center, radius=radius)

        # Rotate the atoms of twisted layer
        orb_pos = layer_twisted.get_positions()
        orb_pos = rotate_coord(orb_pos, angle=angle, center=center)
        orb_pos += np.array([0, 0, shift])
        layer_twisted.positions = orb_pos

    # Merge layers
    quasi_crystal = layer_fixed + layer_twisted

    # Visualization
    vasp_args = {"direct": True, "sort": True}
    write("POSCAR.qs", quasi_crystal, **vasp_args)


def main():
    prim_cell = read("POSCAR", index=-1)
    dim = (33, 33, 1)
    angle = 30 / 180 * math.pi
    center = np.array((2./3, 2./3, 0))
    make_quasi_crystal(prim_cell, dim, angle, center)


if __name__ == "__main__":
    main()
