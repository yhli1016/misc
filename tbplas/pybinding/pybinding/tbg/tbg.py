#! /usr/bin/env python
import math
import time
from typing import List

import numpy as np
from numpy.linalg import norm
from scipy.spatial import cKDTree

import pybinding as pb

import lattice
from utils import Timer


FAIR_PLAY = True


def calc_twist_angle(i: int) -> float:
    """
    Calculate twisting angle according to ref. [1].

    :param i: parameter controlling the twisting angle
    :return: twisting angle in RADIANs, NOT degrees
    """
    cos_ang = (3 * i ** 2 + 3 * i + 0.5) / (3 * i ** 2 + 3 * i + 1)
    return math.acos(cos_ang)


def calc_hetero_lattice(i: int, lat_vec: np.ndarray) -> np.ndarray:
    """
    Calculate Cartesian coordinates of lattice vectors of hetero-structure
    according to ref. [1].

    :param i: parameter controlling the twisting angle
    :param lat_vec: lattice vector
    :return: (3, 3) float64 array, Cartesian coordinates of hetero-structure
        lattice vectors in NANOMETER
    """
    hetero_lattice = np.array([[i, i + 1, 0],
                               [-(i + 1), 2 * i + 1, 0],
                               [0, 0, 1]])
    hetero_lattice = lattice.frac2cart(lat_vec, hetero_lattice)
    return hetero_lattice


def calc_hop(rij: np.ndarray) -> np.ndarray:
    """
    Calculate hopping parameter according to Slater-Koster relation.
    See ref. [2] for the formulae.

    :param rij: (num_rij, 3) array, displacement vector between two orbitals in NM
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
    # Vector operations from now on
    dr = norm(rij, axis=1)
    n = rij[:,2] / dr
    v_pp_pi = - gamma0 * np.exp(q_pi * (1 - dr / a0))
    v_pp_sigma = gamma1 * np.exp(q_sigma * (1 - dr / a1))
    fc = 1 / (1 + np.exp((dr - r_c) / l_c))
    hop = (n**2 * v_pp_sigma + (1 - n**2) * v_pp_pi) * fc
    return hop


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


def make_site_generators(model: pb.Model):
    """
    Make the site generators of orbitals in twisted layer.

    :param model: model for extracting orbital positions, typically has hopping terms
    :return: site generators
    """
    subs = model.system.sublattices
    idx_a = np.flatnonzero(subs == model.lattice.sublattices["A"].alias_id)
    idx_b = np.flatnonzero(subs == model.lattice.sublattices["B"].alias_id)
    pos_a = model.system[idx_a].positions
    pos_b = model.system[idx_b].positions

    on_site = 0.1 if FAIR_PLAY else 0.0
    @pb.site_generator(name="A2", energy=-on_site)
    def add_a():
        return pos_a

    @pb.site_generator(name="B2", energy=on_site)
    def add_b():
        return pos_b

    return add_a, add_b


def make_hop_generators(max_distance: float = 0.75) -> List[pb.hopping_generator]:
    """
    Make intra-layer and inter-layer hopping terms generators.

    :param max_distance: cutoff distance for hopping terms in nm
    :return: list of generators
    """
    d_min, d_max = 0.01, max_distance

    def _kdtree(_pos, _layer1, _layer2):
        kdtree1 = cKDTree(_pos[_layer1])
        kdtree2 = cKDTree(_pos[_layer2])
        coo = kdtree1.sparse_distance_matrix(kdtree2, d_max, output_type='coo_matrix')
        return coo

    @pb.hopping_generator('intra1', energy=-2.7)
    def intra1_gen(x, y, z):
        positions = np.stack([x, y, z], axis=1)
        layer1 = (z == 0)
        coo = _kdtree(positions, layer1, layer1)
        idx = coo.data > d_min
        abs_idx = np.flatnonzero(layer1)
        row, col = abs_idx[coo.row[idx]], abs_idx[coo.col[idx]]
        return row, col

    @pb.hopping_generator('intra2', energy=-2.7)
    def intra2_gen(x, y, z):
        positions = np.stack([x, y, z], axis=1)
        layer2 = (z != 0)
        coo = _kdtree(positions, layer2, layer2)
        idx = coo.data > d_min
        abs_idx = np.flatnonzero(layer2)
        row, col = abs_idx[coo.row[idx]], abs_idx[coo.col[idx]]
        return row, col

    @pb.hopping_generator('inter', energy=-0.1)
    def inter_gen(x, y, z):
        positions = np.stack([x, y, z], axis=1)
        layer1 = (z == 0)
        layer2 = (z != 0)
        coo = _kdtree(positions, layer1, layer2)
        idx = coo.data > d_min
        abs_idx1 = np.flatnonzero(layer1)
        abs_idx2 = np.flatnonzero(layer2)
        row, col = abs_idx1[coo.row[idx]], abs_idx2[coo.col[idx]]
        return row, col

    return intra1_gen, intra2_gen, inter_gen


def make_hop_modifiers() -> List[pb.hopping_energy_modifier]:
    """
    Make modifiers to update hopping terms using Slater-Koster relation.
    """
    @pb.hopping_energy_modifier
    def intra_mod(energy, x1, y1, z1, x2, y2, z2, hop_id):
        rij = np.stack([x2-x1, y2-y1, z2-z1], axis=1)
        intra = (hop_id == 'intra1' or hop_id == "intra2")
        energy[intra] = calc_hop(rij[intra])
        return energy

    @pb.hopping_energy_modifier
    def inter_mod(energy, x1, y1, z1, x2, y2, z2, hop_id):
        rij = np.stack([x2-x1, y2-y1, z2-z1], axis=1)
        inter = (hop_id == 'inter')
        energy[inter] = calc_hop(rij[inter])
        return energy

    return intra_mod, inter_mod


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
    i = 100
    angle = calc_twist_angle(i)
    lat_vec_rot = lattice.rotate_coord(lat_vec, angle)
    orb_pos_rot = lattice.rotate_coord(orb_pos, angle)
    orb_pos_rot[:, 2] += 0.3349

    # Create lattices of two layers
    fixed_lattice = make_pb_lattice(lat_vec, orb_pos)
    twisted_lattice = make_pb_lattice(lat_vec_rot, orb_pos_rot)

    # Make model of top-layer for making modifiers
    hetero_lattice = calc_hetero_lattice(i, lat_vec)
    twisted_model = make_pb_model(hetero_lattice, twisted_lattice)

    tbg = make_pb_model(hetero_lattice,
                        fixed_lattice,
                        make_site_generators(twisted_model),
                        make_hop_generators(),
                        make_hop_modifiers())

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
