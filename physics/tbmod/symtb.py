from collections import defaultdict, namedtuple
from typing import Tuple, Union, List
import math

import sympy as sp
import numpy as np
from scipy.spatial import KDTree


# Type aliases
f_type = Union[int, float, sp.Basic]
c_type = Union[int, float, complex, sp.Basic]
rn_type = Tuple[int, int, int]
pos_type = Tuple[f_type, f_type, f_type]
Orbital = namedtuple("Orbital", ("position", "energy"))
HopTerm = namedtuple("HopTerm", ("rn", "pair", "rij", "distance"))


def cart2frac(lattice_vectors: np.ndarray,
              cartesian_coordinates: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert Cartesian coordinates to fractional coordinates.

    :param lattice_vectors: (3, 3) float64 array
        Cartesian coordinates of lattice vectors
    :param cartesian_coordinates: (num_coord, 3) float64 array
        Cartesian coordinates to convert
    :param origin: float64 array of length 3
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        fractional coordinates in basis of lattice vectors
    """
    if not isinstance(lattice_vectors, np.ndarray):
        lattice_vectors = np.array(lattice_vectors)
    if not isinstance(cartesian_coordinates, np.ndarray):
        cartesian_coordinates = np.array(cartesian_coordinates)
    fractional_coordinates = np.zeros(cartesian_coordinates.shape)
    conversion_matrix = np.linalg.inv(lattice_vectors.T)
    for i, row in enumerate(cartesian_coordinates):
        fractional_coordinates[i] = np.matmul(conversion_matrix,
                                              (row - origin).T)
    return fractional_coordinates


def frac2cart(lattice_vectors: np.ndarray,
              fractional_coordinates: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert fractional coordinates to Cartesian coordinates.

    :param lattice_vectors: (3, 3) float64 array
        Cartesian coordinates of lattice vectors
    :param fractional_coordinates: (num_coord, 3) float64 array
        fractional coordinates to convert in basis of lattice vectors
    :param origin: float64 array of length 3
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        Cartesian coordinates converted from fractional coordinates
    """
    if not isinstance(lattice_vectors, np.ndarray):
        lattice_vectors = np.array(lattice_vectors)
    if not isinstance(fractional_coordinates, np.ndarray):
        fractional_coordinates = np.ndarray(fractional_coordinates)
    cartesian_coordinates = np.zeros(fractional_coordinates.shape)
    conversion_matrix = lattice_vectors.T

    for i, row in enumerate(fractional_coordinates):
        cartesian_coordinates[i] = np.matmul(conversion_matrix, row.T) + origin
    return cartesian_coordinates


def check_conj(hop_ind: Tuple[int, ...], i: int = 0) -> bool:
    """
    Check whether to take the conjugate part of the hopping term.

    :param hop_ind: (r_a, r_b, r_c, orb_i, orb_j), hopping index
    :param i: component index
    :return: whether to take conjugate
    """
    if hop_ind[i] > 0:
        return False
    elif hop_ind[i] < 0:
        return True
    else:
        if i < 2:
            return check_conj(hop_ind, i+1)
        else:
            return hop_ind[3] > hop_ind[4]


class Model:
    """
    Class representing a tight-binding model.

    Attributes
    ----------
    _lattice: (3, 3) float64 np.ndarray
        Cartesian coordinates of lattice vectors
    _orbitals: List[Orbital]
        list of orbital positions and on-site energies
    _hoppings: Dict[Tuple[int, int, int, int, int], c_type]
        keys: cell indices + orbital pairs
        values: hopping energies
    """
    def __init__(self, lattice: np.ndarray = np.eye(3)) -> None:
        """
        :param lattice: (3, 3) float64 np.ndarray
            Cartesian coordinates of lattice vectors
        """
        self._lattice = lattice
        self._orbitals = []
        self._hoppings = dict()

    def add_orbital(self, position: pos_type, energy: f_type = 0) -> None:
        """
        Add a new orbital to the primitive cell.

        :param position: FRACTIONAL coordinate of the orbital
        :param energy: on-site energy of the orbital in eV
        :return: None
        """
        self._orbitals.append(Orbital(position, energy))

    def add_hopping(self, rn: rn_type,
                    orb_i: int,
                    orb_j: int,
                    energy: c_type = 0) -> None:
        """
        Add a new hopping term to the primitive cell, or update an existing
        hopping term.

        :param rn: cell index of the hopping term, i.e. R
        :param orb_i: index of orbital i in <i,0|H|j,R>
        :param orb_j: index of orbital j in <i,0|H|j,R>
        :param energy: hopping integral in eV
        :return: None
        :raises IndexError: if orb_i or orb_j falls out of range
        :raises ValueError: if rn == (0, 0, 0) and orb_i == orb_j
        """
        # Check arguments
        num_orb = len(self._orbitals)
        if not (0 <= orb_i < num_orb):
            raise IndexError(f"orb_i {orb_i} out of range {num_orb}")
        if not (0 <= orb_j < num_orb):
            raise IndexError(f"orb_i {orb_i} out of range {num_orb}")
        if rn == (0, 0, 0) and orb_i == orb_j:
            raise ValueError(f"{rn} {orb_i, orb_j} is an on-site term")

        # Add term
        pair = (orb_i, orb_j)
        if check_conj(rn + pair):
            rn = (-rn[0], -rn[1], -rn[2])
            pair = (orb_j, orb_i)
            energy = energy.conjugate()
        self._hoppings[rn + pair] = energy

    def find_neighbors(self, a_max: int = 0,
                       b_max: int = 0,
                       c_max: int = 0,
                       max_distance: float = 1.0,
                       with_conj: bool = False) -> List[HopTerm]:
        """
        Find neighbours between (0, 0, 0) and nearby cells up to given cutoff
        distance.

        NOTE: only neighbours with distance > 0 will be returned.

        The searching range of nearby cells is:
        [-a_max, a_max] * [-b_max, b_max] * [-c_max, c_max].

        :param a_max: upper bound of range on a-axis
        :param b_max: upper bound of range on b-axis
        :param c_max: upper bound of range on c-axis
        :param max_distance: cutoff distance for hopping terms in the same unit
            as lattice vectors
        :param with_conj: whether to include conjugate terms in the results
        :return: list of neighbors as named tuples
        """
        # Get Cartesian coordinates of orbitals
        pos_r0 = np.array([orb.position for orb in self._orbitals])
        pos_r0 = frac2cart(self._lattice, pos_r0)

        # Prepare for the loop
        tree_r0 = KDTree(pos_r0)
        neighbor_rn = [(ia, ib, ic)
                       for ia in range(-a_max, a_max + 1)
                       for ib in range(-b_max, b_max + 1)
                       for ic in range(-c_max, c_max + 1)]

        # Loop over neighboring cells to search for orbital pairs
        neighbors = []
        for rn in neighbor_rn:
            pos_rn = pos_r0 + np.matmul(rn, self._lattice)
            tree_rn = KDTree(pos_rn)
            dist_matrix = tree_r0.sparse_distance_matrix(tree_rn,
                                                         max_distance=max_distance)
            for pair, distance in dist_matrix.items():
                if distance > 0.0 and (not check_conj(rn + pair) or with_conj):
                    i, j = pair
                    rij = pos_rn[j] - pos_r0[i]
                    neighbors.append(HopTerm(rn, pair, rij, distance))
        neighbors = sorted(neighbors, key=lambda x: x.distance)
        return neighbors

    def print_hk(self, convention: int = 1, with_tril: bool = False) -> None:
        """
        Print analytical Hamiltonian as the function of k-point.

        :param convention: convention for setting up the Hamiltonian
        :param with_tril: whether to print lower triangular hopping terms,
            otherwise print upper triangular hopping terms only
        :return: None
        """
        # Collect on-site terms
        hk = defaultdict(int)
        for i, orb in enumerate(self._orbitals):
            hk[(i, i)] = orb.energy

        # Collect hopping terms
        kpt = [sp.Symbol(_, real=True) for _ in ("ka", "kb", "kc")]
        for hop, energy in self._hoppings.items():
            rn, pair = hop[:3], hop[3:5]
            if convention == 1:
                orb_i, orb_j = pair
                pos_i = self._orbitals[orb_i].position
                pos_j = self._orbitals[orb_j].position
                dr = [rn[_] + pos_j[_] - pos_i[_] for _ in range(3)]
            else:
                dr = rn
            k_dot_r = kpt[0] * dr[0] + kpt[1] * dr[1] + kpt[2] * dr[2]
            phase = 2 * sp.pi * k_dot_r
            hk[pair] += energy * sp.exp(sp.I * phase)
            hk[(pair[1], pair[0])] = hk[pair].conjugate()

        # Print
        for pair, formula in hk.items():
            if formula != 0 and (pair[0] <= pair[1] or with_tril):
                ham_ij = f"ham[{pair[0]}, {pair[1]}]"
                formula = sp.sympify(formula)
                print(f"{ham_ij} = {formula}")


def main():
    lattice = 2.46 * np.array([
        [1.0, 0.0, 0.0],
        [0.5, math.sqrt(3)/2, 0.0],
        [0.0, 0.0, 1.0]
    ])
    f = sp.Rational(1, 3)
    t = sp.Symbol("t", real=True)

    model = Model(lattice)
    model.add_orbital((f, f, 0))
    model.add_orbital((2*f, 2*f, 0))
    model.add_hopping((0, 0, 0), 0, 1, t)
    model.add_hopping((1, 0, 0), 1, 0, t)
    model.add_hopping((0, 1, 0), 1, 0, t)
    model.print_hk(1)
    model.print_hk(2)

    neighbors = model.find_neighbors(a_max=2, b_max=2, max_distance=2.0)
    for term in neighbors:
        print(term.rn, term.pair, term.distance)


if __name__ == "__main__":
    main()
