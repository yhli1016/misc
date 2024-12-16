from collections import namedtuple
from typing import Tuple, List, Union

import sympy as sp
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

from .utils import f_type, c_type, pos_type, rn_type, frac2cart
from .visual import ModelViewer


__all__ = ["Model"]


# Type aliases
Orbital = namedtuple("Orbital", ("position", "energy"))
HopTerm = namedtuple("HopTerm", ("rn", "pair", "rij", "distance"))


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
    _lattice: Union[sp.Matrix, np.ndarray]
        Cartesian coordinates of lattice vectors
    _orbitals: List[Orbital]
        list of orbital positions and on-site energies
    _hoppings: Dict[Tuple[int, int, int, int, int], c_type]
        keys: cell indices + orbital pairs
        values: hopping energies
    """
    def __init__(self, lattice: Union[sp.Matrix, np.ndarray] = None) -> None:
        """
        :param lattice: (3, 3) sympy or numpy array
            Cartesian coordinates of lattice vectors
        """
        if lattice is None:
            self._lattice = sp.eye(3)
        else:
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
        num_orb = self.num_orb
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

    def get_hk(self, convention: int = 1) -> sp.Matrix:
        """
        Get analytical Hamiltonian matrix as the function of k-point.

        :param convention: convention for setting up the Hamiltonian
        :return: analytical Hamiltonian matrix
        """
        # Collect on-site terms
        hk = sp.zeros(self.num_orb)
        for i, orb in enumerate(self._orbitals):
            hk[i, i] = orb.energy

        # Collect hopping terms
        kpt = [sp.Symbol(_, real=True) for _ in ("ka", "kb", "kc")]
        for hop, energy in self._hoppings.items():
            rn, orb_i, orb_j = hop[:3], hop[3], hop[4]
            if convention == 1:
                pos_i = self._orbitals[orb_i].position
                pos_j = self._orbitals[orb_j].position
                dr = [rn[_] + pos_j[_] - pos_i[_] for _ in range(3)]
            else:
                dr = rn
            k_dot_r = kpt[0] * dr[0] + kpt[1] * dr[1] + kpt[2] * dr[2]
            phase = 2 * sp.pi * k_dot_r
            hij = energy * sp.exp(sp.I * phase)
            # DO NOT use hk[i, j] = hk[j, i].conjugate(), otherwise the diagonal
            # terms will be wrong!
            hk[orb_i, orb_j] += hij
            hk[orb_j, orb_i] += hij.conjugate()
        return hk

    def print_hk(self, convention: int = 1) -> None:
        """
        Print analytical hamiltonian

        :param convention: convention for setting up the Hamiltonian
        :return:
        """
        gauge = "Atomic" if convention == 1 else "Lattice"
        print(f"{gauge} gauge:")
        hk = self.get_hk(convention)
        for ii in range(self.num_orb):
            for jj in range(self.num_orb):
                print(f"ham[{ii}, {jj}] = {hk[ii, jj]}")

    def print_cxx(self) -> None:
        """
        Print c++ code for constructing model.
        :return: None
        """
        print("// Units assumed to be ANGSTROM/eV.\n")

        # Lattice vectors and origin
        print("// Lattice vectors and origin.")
        print("Eigen::Matrix3d lat_vec{")
        for i in range(3):
            print("{", end="")
            for j in range(2):
                print(" ", self._lattice[i, j], ",", end="")
            if i < 2:
                print(" ", self._lattice[i, j], "},")
            else:
                print(" ", self._lattice[i, j], "}")
        print("};")
        print("Eigen::Vector3d origin(0.0, 0.0, 0.0);\n")

        # Create the primitive cell and set orbitals
        print("// Create the primitive cell and set orbitals.")
        print(f"PrimitiveCell<complex_t> prim_cell({self.num_orb}, lat_vec,"
              f" origin, tbplas::base::ANG);")
        for i, orbital in enumerate(self._orbitals):
            pos = orbital.position
            eng = orbital.energy
            print(f"prim_cell.set_orbital{i, pos[0], pos[1], pos[2], eng, 0};")
        print()

        # Add hopping terms
        print("// Add hopping terms.")
        for rn, eng in self._hoppings.items():
            if eng != 0:
                print(f"prim_cell.add_hopping{rn[0], rn[1], rn[2], rn[3], rn[4], eng};")

    def plot(self, fig_name: str = None,
             fig_size: Tuple[float, float] = None,
             fig_dpi: int = 300,
             with_orbitals: bool = True,
             with_cells: bool = True,
             with_conj: bool = True,
             orb_color: List[str] = None,
             hop_as_arrows: bool = True,
             hop_color: str = "r",
             view: str = "ab") -> None:
        """
        Plot lattice vectors, orbitals, and hopping terms.

        If figure name is given, save the figure to file. Otherwise, show it on
        the screen.

        :param fig_name: file name to which the figure will be saved
        :param fig_size: width and height of the figure
        :param fig_dpi: resolution of the figure file
        :param with_orbitals: whether to plot orbitals as filled circles
        :param with_cells: whether to plot borders of primitive cells
        :param with_conj: whether to plot conjugate hopping terms as well
        :param orb_color: colors of the orbitals
        :param hop_as_arrows: whether to plot hopping terms as arrows
            If true, hopping terms will be plotted as arrows using axes.arrow()
            method. Otherwise, they will be plotted as lines using
            LineCollection. The former is more intuitive but much slower.
        :param hop_color: color of hopping terms
        :param view: kind of view point
        :returns: None
        :raises ValueError: if view is illegal
        """
        # Assemble arrays
        hop_ind = np.array([_ for _ in self._hoppings.keys()])
        orb_pos = np.array([orb.position for orb in self._orbitals])
        orb_pos = frac2cart(self._lattice, orb_pos)
        origin = np.zeros(3)
        dr = np.zeros((self.num_hop, 3), dtype=np.float64)
        for i_h, ind in enumerate(hop_ind):
            orb_i, orb_j = ind.item(3), ind.item(4)
            rn = np.matmul(ind[0:3], self._lattice)
            dr[i_h] = orb_pos[orb_j] + rn - orb_pos[orb_i]

        # Initialize visualizer
        fig, axes = plt.subplots(figsize=fig_size)
        axes.set_aspect('equal')
        viewer = ModelViewer(axes, self._lattice, origin, view)

        # Determine the range of rn
        rn_range = np.zeros((3, 2), dtype=np.int32)
        if self.num_hop > 0:
            for i in range(3):
                ri_min = hop_ind[:, i].min()
                ri_max = hop_ind[:, i].max()
                if with_conj:
                    rn_range[i, 0] = min([ri_min, ri_max, -ri_min, -ri_max])
                    rn_range[i, 1] = max([ri_min, ri_max, -ri_min, -ri_max])
                else:
                    rn_range[i, 0] = ri_min
                    rn_range[i, 1] = ri_max
        ra_min, ra_max = rn_range.item(0, 0), rn_range.item(0, 1)
        rb_min, rb_max = rn_range.item(1, 0), rn_range.item(1, 1)
        rc_min, rc_max = rn_range.item(2, 0), rn_range.item(2, 1)

        # Plot orbitals
        if orb_color is None:
            orb_color = ['b' for _ in range(self.num_orb)]
        if self.num_orb > 0:
            if with_orbitals:
                for i_a in range(ra_min, ra_max+1):
                    for i_b in range(rb_min, rb_max+1):
                        for i_c in range(rc_min, rc_max+1):
                            center = np.matmul((i_a, i_b, i_c), self._lattice)
                            pos_rn = orb_pos + center
                            viewer.scatter(pos_rn, s=100, c=orb_color)

        # Plot hopping terms
        if self.num_hop > 0:
            hop_i = hop_ind[:, 3]
            hop_j = hop_ind[:, 4]
            arrow_args = {"color": hop_color, "length_includes_head": True,
                          "width": 0.02, "head_width": 0.2, "fill": False}
            for i_h in range(hop_i.shape[0]):
                # Original term
                pos_i = orb_pos[hop_i.item(i_h)]
                pos_j = pos_i + dr[i_h]
                if hop_as_arrows:
                    viewer.plot_arrow(pos_i, pos_j, **arrow_args)
                else:
                    viewer.add_line(pos_i, pos_j)

                # Conjugate term
                if with_conj:
                    pos_j = orb_pos[hop_j.item(i_h)]
                    pos_i = pos_j - dr[i_h]
                    if hop_as_arrows:
                        viewer.plot_arrow(pos_j, pos_i, **arrow_args)
                    else:
                        viewer.add_line(pos_j, pos_i)
            if not hop_as_arrows:
                viewer.plot_line(color=hop_color)

        # Plot cells
        if with_cells:
            if view in ("ab", "ba"):
                viewer.add_grid(ra_min, ra_max + 1, rb_min, rb_max + 1)
            elif view in ("bc", "cb"):
                viewer.add_grid(rb_min, rb_max + 1, rc_min, rc_max + 1)
            else:
                viewer.add_grid(ra_min, ra_max + 1, rc_min, rc_max + 1)
            viewer.plot_grid(color="k", linestyle=":")
            viewer.plot_lat_vec(color="k", length_includes_head=True,
                                width=0.05, head_width=0.2)

        # Hide spines and ticks.
        for key in ("top", "bottom", "left", "right"):
            axes.spines[key].set_visible(False)
        axes.set_xticks([])
        axes.set_yticks([])
        fig.tight_layout()
        plt.autoscale()
        if fig_name is not None:
            plt.savefig(fig_name, dpi=fig_dpi)
        else:
            plt.show()
        plt.close()

    @property
    def num_orb(self) -> int:
        """Get the number of orbitals in the model."""
        return len(self._orbitals)

    @property
    def num_hop(self) -> int:
        """Get the number of hopping terms in the model."""
        return len(self._hoppings)
