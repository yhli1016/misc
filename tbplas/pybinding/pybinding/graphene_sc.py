#! /usr/bin/env python
from typing import List

import numpy as np


import pybinding as pb

from utils import Timer


FAIR_PLAY = True


def make_pbc(a1: bool, a2: bool, a3: bool) -> pb.translational_symmetry:
    return pb.translational_symmetry(a1=a1, a2=a2, a3=a3)


def make_pb_lattice(lat_vec: np.ndarray, orb_pos: np.ndarray) -> pb.Lattice:
    """
    Make an lattice with no hopping terms, which will be generated with hop generators.

    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors in nm, with each row being a vector
    :param orb_pos: (num_orb, 3) float64 array
        Cartesian coordinates of lattice vectors in nm, with each row being a vector
    :return: lattice with no hopping terms
    """
    lattice = pb.Lattice(a1=lat_vec[0], a2=lat_vec[1], a3=lat_vec[2])
    on_site = 0.1 if FAIR_PLAY else 0.0
    lattice.add_sublattices(
        ('A', orb_pos[0], -on_site),
        ('B', orb_pos[1], on_site))
    lattice.add_hoppings(
        ([0, 0], 'A', 'B', -2.7),
        ([1, 0], 'B', 'A', -2.7),
        ([0, 1], 'B', 'A', -2.7)
    )
    return lattice


def make_pb_model(hetero_lattice: np.ndarray, pb_lattice: pb.Lattice, *args) -> pb.Model:
    """
    Make a tbg model.

    :param hetero_lattice: (3, 3) float64 array
        Cartesian coordinates of hetero-structure lattice vectors in nm
    :param pb_lattice: pb.Lattice to create the model
    :param args: list of pybinding modifiers
    :return: pybinding model
    """
    a1, a2 = hetero_lattice[0], hetero_lattice[1]
    shape = pb.Polygon([[0.0, 0.0, 0.0], a1, a1+a2, a2])
    pbc = pb.translational_symmetry(a1=False, a2=False, a3=True)
    return pb.Model(pb_lattice, shape, pbc, *args)


def main():
    # Cartesian coordinates in nm
    lat_vec = np.array([
        [2.46000000e-01, 0.00000000e+00, 0.00000000e+00],
        [1.23000000e-01, 2.13042249e-01, 0.00000000e+00],
        [6.12323400e-17, 3.53525080e-17, 1.00000000e+00]])
    orb_pos = np.array([
        [0.        , 0.        , 0.        ],
        [0.123     , 0.07101408, 0.        ]])

    # Rotate and shift coordinates
    sc_lat_vec = lat_vec.copy()
    sc_lat_vec[0] *= 1024
    sc_lat_vec[1] *= 1024

    # Create lattices of two layers
    fixed_lattice = make_pb_lattice(lat_vec, orb_pos)

    # Make model of top-layer for making modifiers
    tbg = make_pb_model(sc_lat_vec, fixed_lattice)

    #tbg.plot()
    #plt.show()
    if FAIR_PLAY:
        ham = tbg.hamiltonian
        print(ham.shape[0], ham.nnz)


if __name__ == "__main__":
    timer = Timer()
    for i in range(5):
        timer.tic("tbg")
        main()
        timer.toc("tbg")
        timer.report_time()
