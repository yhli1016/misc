"""Functions and classes for advanced modeling."""

from typing import Tuple, Union
import math

import sympy as sp
import numpy as np


__all__ = ["f_type", "c_type", "rn_type", "pos_type", "cart2frac", "frac2cart",
           "SK", "SOCTable"]


# Type aliases
f_type = Union[int, float, sp.Basic]
c_type = Union[int, float, complex, sp.Basic]
rn_type = Tuple[int, int, int]
pos_type = Tuple[f_type, f_type, f_type]


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


class SK:
    """
    Class for evaluating hopping integrals using Slater-Koster formula.

    The maximum supported angular momentum is l=2 (d orbitals).
    Reference: https://journals.aps.org/pr/abstract/10.1103/PhysRev.94.1498

    The reason why we make orbital labels as attributes is to avoid misspelling,
    which is likely to happen as we have to repeat them many times.
    """
    def __init__(self) -> None:
        self._s = "s"
        self._px = "px"
        self._py = "py"
        self._pz = "pz"
        self._dxy = "dxy"
        self._dyz = "dyz"
        self._dzx = "dzx"
        self._dx2_y2 = "dx2-y2"
        self._dz2 = "dz2"
        self._p_labels = {self._px, self._py, self._pz}
        self._d_labels = {self._dxy, self._dyz, self._dzx, self._dx2_y2,
                          self._dz2}
        self._sqrt3 = math.sqrt(3)
        self._half_sqrt3 = self._sqrt3 * 0.5

    def _check_p_labels(self, *labels: str) -> None:
        """
        Check the sanity of labels of p orbitals.

        :param labels: labels to check
        :return: None
        :raises ValueError: if any label is not in self.p_labels
        """
        for label in labels:
            if label not in self._p_labels:
                raise ValueError(f"Illegal label: {label}")

    def _check_d_labels(self, *labels: str) -> None:
        """
        Check the sanity of labels of d orbitals.

        :param labels: labels to check
        :return: None
        :raises ValueError: if any label is not in self.d_labels
        """
        for label in labels:
            if label not in self._d_labels:
                raise ValueError(f"Illegal label: {label}")

    @staticmethod
    def _perm_vector(vector: np.ndarray, x_new: str) -> np.ndarray:
        """
        Permute a given vector according to the new x_axis.

        :param vector: vector to permute
        :param x_new: label of the new x_axis
        :return: permuted vector
        :raises ValueError: if x_new is not in 'x', 'y', 'z'
        """
        if x_new == "x":
            return vector
        elif x_new == "y":
            return vector[[1, 2, 0]]
        elif x_new == "z":
            return vector[[2, 0, 1]]
        else:
            raise ValueError(f"Illegal x_new {x_new}")

    @staticmethod
    def _remap_label(label: str, x_new: str) -> str:
        """
        Remap orbital label after permutation.

        :param label: orbital label to remap
        :param x_new: index of the new 'x' direction
        :return: remapped orbital label
        :raises ValueError: if x_new is not in 'x', 'y', 'z'
        """
        if x_new == "x":
            return label
        else:
            if x_new == "y":
                # keys: x, y, z, values: x_new, y_new, z_new
                map_table = {"x": "z", "y": "x", "z": "y"}
            elif x_new == "z":
                # keys: x, y, z, values: x_new, y_new, z_new
                map_table = {"x": "y", "y": "z", "z": "x"}
            else:
                raise ValueError(f"Illegal new_x {x_new}")
            new_label = label[0]
            for c in label[1:]:
                new_label += map_table[c]
            return new_label

    @staticmethod
    def _eval_dir_cos(r: np.ndarray) -> Tuple[float, float, float]:
        """
        Evaluate direction cosines for given displacement vector.

        :param r: Cartesian coordinates of the displacement vector
        :return: the direction cosines along x, y, and z directions
        :raises ValueError: if the norm of r is too small
        """
        norm = np.linalg.norm(r)
        if norm <= 1.0e-15:
            raise ValueError("Norm of displacement vector too small")
        dir_cos = r / norm
        l, m, n = dir_cos.item(0), dir_cos.item(1), dir_cos.item(2)
        return l, m, n

    @staticmethod
    def ss(v_sss: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <s,0|H|s,r>.

        :param v_sss: V_ss_sigma
        :return: hopping integral
        """
        return v_sss

    def sp(self, r: np.ndarray,
           label_p: str = "px",
           v_sps: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <s,0|H|p,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_p: label of p orbital
        :param v_sps: V_sp_sigma
        :return: hopping integral
        :raises ValueError: if label_p is not in self.p_labels
        """
        self._check_p_labels(label_p)
        l, m, n = self._eval_dir_cos(r)
        if label_p == self._px:
            t = l * v_sps
        elif label_p == self._py:
            t = m * v_sps
        else:
            t = n * v_sps
        return t

    def sd(self, r: np.ndarray,
           label_d: str = "dxy",
           v_sds: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <s,0|H|d,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_d: label of d orbital
        :param v_sds: V_sd_sigma
        :return: hopping integral
        :raises ValueError: if label_d is not in self.d_labels
        """
        self._check_d_labels(label_d)

        # Permute the coordinates
        if label_d == self._dyz:
            x_new = "y"
        elif label_d == self._dzx:
            x_new = "z"
        else:
            x_new = "x"
        r = self._perm_vector(r, x_new)
        label_d = self._remap_label(label_d, x_new)

        # Evaluate the hopping integral
        l, m, n = self._eval_dir_cos(r)
        if label_d == self._dxy:
            t = self._sqrt3 * l * m * v_sds
        elif label_d == self._dx2_y2:
            t = self._half_sqrt3 * (l ** 2 - m ** 2) * v_sds
        elif label_d == self._dz2:
            t = (n ** 2 - 0.5 * (l ** 2 + m ** 2)) * v_sds
        else:
            raise ValueError(f"Undefined label pair s {label_d}")
        return t

    def pp(self, r: np.ndarray,
           label_i: str = "px",
           label_j: str = "px",
           v_pps: c_type = 0,
           v_ppp: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <p_i,0|H|p_j,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_i: label of p_i orbital
        :param label_j: label of p_j orbital
        :param v_pps: V_pp_sigma
        :param v_ppp: V_pp_pi
        :return: hopping integral
        :raises ValueError: if label_i or label_j is not in self.p_labels
        """
        self._check_p_labels(label_i, label_j)

        # Permute the coordinates
        x_new = label_i[1]
        r = self._perm_vector(r, x_new)
        label_i = self._remap_label(label_i, x_new)
        label_j = self._remap_label(label_j, x_new)

        # After permutation, label_i will always be px.
        # Otherwise, something must be wrong.
        if label_i != self._px:
            raise ValueError(f"Undefined label pair {label_i} {label_j}")

        # The minimal hopping table in the reference.
        l, m, n = self._eval_dir_cos(r)
        if label_j == self._px:
            t = l ** 2 * v_pps + (1 - l ** 2) * v_ppp
        elif label_j == self._py:
            t = l * m * (v_pps - v_ppp)
        else:
            t = l * n * (v_pps - v_ppp)
        return t

    def pd(self, r: np.ndarray,
           label_p: str = "px",
           label_d: str = "dxy",
           v_pds: c_type = 0,
           v_pdp: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <p,0|H|d,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_p: label of p orbital
        :param label_d: label of d orbital
        :param v_pds: V_pd_sigma
        :param v_pdp: V_pd_pi
        :return: hopping integral
        :raises ValueError: if label_p is not in self.p_labels or
            label_d is not in self.d_labels
        """
        self._check_p_labels(label_p)
        self._check_d_labels(label_d)

        # Permute coordinates
        perm_labels = (self._dxy, self._dyz, self._dzx)
        if label_d in perm_labels:
            x_new = label_p[1]
            r = self._perm_vector(r, x_new)
            label_p = self._remap_label(label_p, x_new)
            label_d = self._remap_label(label_d, x_new)

        # Evaluate hopping integral
        l, m, n = self._eval_dir_cos(r)
        l2, m2, n2 = l ** 2, m ** 2, n ** 2
        l2_p_m2 = l2 + m2
        l2_m_m2 = l2 - m2
        sqrt3 = self._sqrt3
        sqrt3_2 = self._half_sqrt3

        if label_p == self._px:
            if label_d == self._dxy:
                t = sqrt3 * l2 * m * v_pds + m * (1 - 2 * l2) * v_pdp
            elif label_d == self._dyz:
                t = l * m * n * (sqrt3 * v_pds - 2 * v_pdp)
            elif label_d == self._dzx:
                t = sqrt3 * l2 * n * v_pds + n * (1 - 2 * l2) * v_pdp
            elif label_d == self._dx2_y2:
                t = sqrt3_2 * l * l2_m_m2 * v_pds + l * (1 - l2_m_m2) * v_pdp
            else:
                t = l * (n2 - 0.5 * l2_p_m2) * v_pds - sqrt3 * l * n2 * v_pdp
        elif label_p == self._py:
            if label_d == self._dx2_y2:
                t = sqrt3_2 * m * l2_m_m2 * v_pds - m * (1 + l2_m_m2) * v_pdp
            elif label_d == self._dz2:
                t = m * (n2 - 0.5 * l2_p_m2) * v_pds - sqrt3 * m * n2 * v_pdp
            else:
                raise ValueError(f"Undefined label pair {label_p} {label_d}")
        else:
            if label_d == self._dx2_y2:
                t = sqrt3_2 * n * l2_m_m2 * v_pds - n * l2_m_m2 * v_pdp
            elif label_d == self._dz2:
                t = n * (n2 - 0.5 * l2_p_m2) * v_pds + sqrt3 * n * l2_p_m2 * v_pdp
            else:
                raise ValueError(f"Undefined label pair {label_p} {label_d}")
        return t

    def dd(self, r: np.ndarray,
           label_i: str = "dxy",
           label_j: str = "dxy",
           v_dds: c_type = 0,
           v_ddp: c_type = 0,
           v_ddd: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <d_i,0|H|d_j,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_i: label of d_i orbital
        :param label_j: label of d_j orbital
        :param v_dds: V_dd_sigma
        :param v_ddp: V_dd_pi
        :param v_ddd: V_dd_delta
        :return: hopping integral
        :raises ValueError: if label_i or label_j is not in self.d_labels
        """
        self._check_d_labels(label_i, label_j)

        # Number the orbitals such that we can filter the diagonal terms
        # The order of the orbitals strictly follows the reference.
        # DO NOT CHANGE IT UNLESS YOU KNOW WHAT YOU ARE DOING!
        d_labels = (self._dxy, self._dyz, self._dzx, self._dx2_y2, self._dz2)
        id_i = d_labels.index(label_i)
        id_j = d_labels.index(label_j)

        if id_i > id_j:
            t = self.dd(r=-r, label_i=label_j, label_j=label_i, v_dds=v_dds,
                        v_ddp=v_ddp, v_ddd=v_ddd).conjugate()
        else:
            # Permute the coordinates if essential
            if label_i == self._dyz and label_j in (self._dyz, self._dzx):
                x_new = "y"
            elif label_i == self._dzx and label_j == self._dzx:
                x_new = "z"
            else:
                x_new = "x"
            r = self._perm_vector(r, x_new)
            label_i = self._remap_label(label_i, x_new)
            label_j = self._remap_label(label_j, x_new)

            # Evaluate hopping integral
            l, m, n = self._eval_dir_cos(r)
            l2, m2, n2 = l ** 2, m ** 2, n ** 2
            l2_p_m2 = l2 + m2
            l2_m_m2 = l2 - m2
            lm, mn, nl = l * m, m * n, n * l
            l2m2 = l2 * m2
            sqrt3 = self._sqrt3

            if label_i == self._dxy:
                if label_j == self._dxy:
                    factor = 1.0
                    t1 = 3 * l2m2
                    t2 = l2_p_m2 - 4 * l2m2
                    t3 = n2 + l2m2
                elif label_j == self._dyz:
                    factor = nl
                    t1 = 3 * m2
                    t2 = 1 - 4 * m2
                    t3 = m2 - 1
                elif label_j == self._dzx:
                    factor = mn
                    t1 = 3 * l2
                    t2 = 1 - 4 * l2
                    t3 = l2 - 1
                elif label_j == self._dx2_y2:
                    factor = lm * l2_m_m2
                    t1 = 1.5
                    t2 = -2
                    t3 = 0.5
                else:
                    factor = sqrt3 * lm
                    t1 = n2 - 0.5 * l2_p_m2
                    t2 = -2 * n2
                    t3 = 0.5 * (1 + n2)
            elif label_i == self._dyz:
                if label_j == self._dx2_y2:
                    factor = mn
                    t1 = 1.5 * l2_m_m2
                    t2 = -(1 + 2 * l2_m_m2)
                    t3 = 1 + 0.5 * l2_m_m2
                elif label_j == self._dz2:
                    factor = sqrt3 * mn
                    t1 = n2 - 0.5 * l2_p_m2
                    t2 = l2_p_m2 - n2
                    t3 = -0.5 * l2_p_m2
                else:
                    raise ValueError(f"Undefined label pair {label_i} {label_j}")
            elif label_i == self._dzx:
                if label_j == self._dx2_y2:
                    factor = nl
                    t1 = 1.5 * l2_m_m2
                    t2 = 1 - 2 * l2_m_m2
                    t3 = -1 * (1 - 0.5 * l2_m_m2)
                elif label_j == self._dz2:
                    factor = sqrt3 * nl
                    t1 = n2 - 0.5 * l2_p_m2
                    t2 = l2_p_m2 - n2
                    t3 = -0.5 * l2_p_m2
                else:
                    raise ValueError(f"Undefined label pair {label_i} {label_j}")
            elif label_i == self._dx2_y2:
                if label_j == self._dx2_y2:
                    factor = 1
                    t1 = 0.75 * l2_m_m2 ** 2
                    t2 = l2_p_m2 - l2_m_m2 ** 2
                    t3 = n2 + 0.25 * l2_m_m2 ** 2
                elif label_j == self._dz2:
                    factor = sqrt3 * l2_m_m2
                    t1 = 0.5 * (n2 - 0.5 * l2_p_m2)
                    t2 = -n2
                    t3 = 0.25 * (1 + n2)
                else:
                    raise ValueError(f"Undefined label pair {label_i} {label_j}")
            else:
                if label_j == self._dz2:
                    factor = 1
                    t1 = (n2 - 0.5 * l2_p_m2) ** 2
                    t2 = 3 * n2 * l2_p_m2
                    t3 = 0.75 * l2_p_m2 ** 2
                else:
                    raise ValueError(f"Undefined label pair {label_i} {label_j}")
            t = factor * (t1 * v_dds + t2 * v_ddp + t3 * v_ddd)
        return t

    def ps(self, r: np.ndarray,
           label_p: str = "px",
           v_sps: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <p,0|H|s,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_p: label of p orbital
        :param v_sps: V_sp_sigma
        :return: hopping integral
        :raises ValueError: if label_p is not in self.p_labels
        """
        return self.sp(r=-r, label_p=label_p, v_sps=v_sps).conjugate()

    def ds(self, r: np.ndarray,
           label_d: str = "dxy",
           v_sds: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <d,0|H|s,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_d: label of d orbital
        :param v_sds: V_sd_sigma
        :return: hopping integral
        :raises ValueError: if label_d is not in self.d_labels
        """
        return self.sd(r=-r, label_d=label_d, v_sds=v_sds).conjugate()

    def dp(self, r: np.ndarray,
           label_p: str = "px",
           label_d: str = "dxy",
           v_pds: c_type = 0,
           v_pdp: c_type = 0) -> c_type:
        """
        Evaluate the hopping integral <d,0|H|p,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_p: label of p orbital
        :param label_d: label of d orbital
        :param v_pds: V_pd_sigma
        :param v_pdp: V_pd_pi
        :return: hopping integral
        :raises ValueError: if label_p is not in self.p_labels or
            label_d is not in self.d_labels
        """
        return self.pd(r=-r, label_p=label_p, label_d=label_d, v_pds=v_pds,
                       v_pdp=v_pdp).conjugate()

    def eval(self, r: np.ndarray,
             label_i: str = "s",
             label_j: str = "s",
             v_sss: c_type = 0, v_sps: c_type = 0, v_sds: c_type = 0,
             v_pps: c_type = 0, v_ppp: c_type = 0,
             v_pds: c_type = 0, v_pdp: c_type = 0,
             v_dds: c_type = 0, v_ddp: c_type = 0, v_ddd: c_type = 0) -> c_type:
        """
        Common interface for evaluating hopping integral <i,0|H|j,r>.

        :param r: Cartesian coordinates of the displacement vector
        :param label_i: label for orbital i
        :param label_j: label for orbital j
        :param v_sss: V_ss_sigma
        :param v_sps: V_sp_sigma
        :param v_sds: V_sd_sigma
        :param v_pps: V_pp_sigma
        :param v_ppp: V_pp_pi
        :param v_pds: V_pd_sigma
        :param v_pdp: V_pd_pi
        :param v_dds: V_dd_sigma
        :param v_ddp: V_dd_pi
        :param v_ddd: V_dd_delta
        :return: hopping integral
        :raises ValueError: if label_i or label_j is not in predefined labels
        """
        if label_i == self._s:
            if label_j == self._s:
                t = self.ss(v_sss=v_sss)
            elif label_j in self._p_labels:
                t = self.sp(r=r, label_p=label_j, v_sps=v_sps)
            else:
                t = self.sd(r=r, label_d=label_j, v_sds=v_sds)
        elif label_i in self._p_labels:
            if label_j == self._s:
                t = self.ps(r=r, label_p=label_i, v_sps=v_sps)
            elif label_j in self._p_labels:
                t = self.pp(r=r, label_i=label_i, label_j=label_j,
                            v_pps=v_pps, v_ppp=v_ppp)
            else:
                t = self.pd(r=r, label_p=label_i, label_d=label_j,
                            v_pds=v_pds, v_pdp=v_pdp)
        else:
            if label_j == self._s:
                t = self.ds(r=r, label_d=label_i, v_sds=v_sds)
            elif label_j in self._p_labels:
                t = self.dp(r=r, label_d=label_i, label_p=label_j,
                            v_pds=v_pds, v_pdp=v_pdp)
            else:
                t = self.dd(r=r, label_i=label_i, label_j=label_j,
                            v_dds=v_dds, v_ddp=v_ddp, v_ddd=v_ddd)
        return t


class SOCTable:
    """
    Hard-coded spin-orbital coupling term table, taking from
    https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.245412
    adapted from the 'SOCTable2' class of TBPLaS.

    Attributes
    ----------
    _data: Dict[Tuple[str, str], sp.Matrix]
        soc coupling matrices
    _orbital_labels: Set[str]
        labels of atomic orbitals
    _spin_labels: Set[str]
        directions of spins
    """
    def __init__(self):
        half = sp.Rational(1, 2)
        s_x = half * sp.Matrix([[0, 1], [1, 0]])
        s_y = half * sp.Matrix([[0, -sp.I], [sp.I, 0]])
        s_z = half * sp.Matrix([[1, 0], [0, -1]])
        self._data = {
            ('px', 'py'): -sp.I * s_z,
            ('px', 'pz'): sp.I * s_y,
            ('py', 'pz'): -sp.I * s_x,
            ('dxy', 'dx2-y2'): 2 * sp.I * s_z,
            ('dxy', 'dzx'): -sp.I * s_x,
            ('dxy', 'dyz'): sp.I * s_y,
            ('dx2-y2', 'dzx'): sp.I * s_y,
            ('dx2-y2', 'dyz'): sp.I * s_x,
            ('dzx', 'dyz'): -sp.I * s_z,
            ('dzx', 'dz2'): sp.I * sp.sqrt(3) * s_y,
            ('dyz', 'dz2'): -sp.I * sp.sqrt(3) * s_x,
        }
        # Restore the remaining terms according to antisymmetric relation
        data2 = dict()
        for key, value in self._data.items():
            data2[(key[1], key[0])] = -value
        self._data.update(data2)
        self._orbital_labels = {"s", "px", "py", "pz",
                                "dxy", "dx2-y2", "dyz", "dzx", "dz2"}
        self._spin_labels = {"up", "down"}

    def eval(self, label_i: str = "s",
             spin_i: str = "up",
             label_j: str = "s",
             spin_j: str = "down") -> c_type:
        """
        Evaluate the matrix element <i,s_i|l*s|j,s_j>.

        :param label_i: orbital label of bra
        :param spin_i: spin direction of bra
        :param label_j: orbital label of ket
        :param spin_j: spin direction of ket
        :return: matrix element in h_bar**2
        :raises ValueError: if orbital labels or spin directions are illegal
        """
        label_idx = (label_i, label_j)
        spin_idx = (spin_i, spin_j)
        for label in label_idx:
            if label not in self._orbital_labels:
                raise ValueError(f"Illegal orbital label {label}")
        for spin in spin_idx:
            if spin not in self._spin_labels:
                raise ValueError(f"Illegal spin direction {spin}")
        try:
            soc_mat = self._data[label_idx]
        except KeyError:
            soc_mat = sp.zeros(2)
        s_i = 0 if spin_i == "up" else 1
        s_j = 0 if spin_j == "up" else 1
        return soc_mat[s_i, s_j]
