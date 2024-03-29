from collections import defaultdict
from typing import Tuple, Union

import sympy as sp


# Type aliases
f_type = Union[int, float, sp.Basic]
c_type = Union[int, float, complex, sp.Basic]
rn_type = Tuple[int, int, int]
pair_type = Tuple[int, int]
pos_type = Tuple[f_type, f_type, f_type]


def invert_rn(rn: rn_type, i: int = 0) -> bool:
    """
    Check if the cell index should be inverted.

    :param rn: (r_a, r_b, r_c), cell index
    :param i: component index
    :return: whether to invert the cell index
    """
    if rn[i] > 0:
        return False
    elif rn[i] < 0:
        return True
    else:
        if i < 2:
            return invert_rn(rn, i+1)
        else:
            return False


def norm_keys(rn: Tuple[int, int, int],
              orb_i: int,
              orb_j: int) -> Tuple[rn_type, pair_type, bool]:
    """
    Normalize cell index and orbital pair into permitted keys.

    For IntraHopping, it should check whether to take the conjugation and
    return the status in conj.

    :param rn: (r_a, r_b, r_c), cell index
    :param orb_i: orbital index or bra
    :param orb_j: orbital index of ket
    :return: (rn, pair, conj)
        where rn is the normalized cell index,
        pair is the normalized orbital pair,
        conj is the flag of whether to take the conjugate of hopping energy
    """
    if invert_rn(rn):
        rn = (-rn[0], -rn[1], -rn[2])
        pair = (orb_j, orb_i)
        conj = True
    elif rn == (0, 0, 0) and orb_i > orb_j:
        rn = rn
        pair = (orb_j, orb_i)
        conj = True
    else:
        rn = rn
        pair = (orb_i, orb_j)
        conj = False
    return rn, pair, conj


class Model:
    """
    Class representing a tight-binding model.

    Attributes
    ----------
    _orbitals: List[Tuple[pos_type, f_type]]
        list of orbital positions and on-site energies
    _hoppings: Dict[Tuple[int, int, int, int, int], c_type]
        keys: cell indices + orbital pairs
        values: hopping energies
    """
    def __init__(self):
        self._orbitals = []
        self._hoppings = dict()

    def add_orbital(self, position: pos_type, energy: f_type = 0) -> None:
        """
        Add a new orbital to the primitive cell.

        :param position: FRACTIONAL coordinate of the orbital
        :param energy: on-site energy of the orbital in eV
        :return: None
        """
        self._orbitals.append((position, energy))

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
        num_orb = len(self._orbitals)
        if not (0 <= orb_i < num_orb):
            raise IndexError(f"orb_i {orb_i} out of range {num_orb}")
        if not (0 <= orb_j < num_orb):
            raise IndexError(f"orb_i {orb_i} out of range {num_orb}")
        if rn == (0, 0, 0) and orb_i == orb_j:
            raise ValueError(f"{rn} {orb_i, orb_j} is an on-site term")
        rn, pair, conj = norm_keys(rn, orb_i, orb_j)
        if conj:
            energy = energy.conjugate()
        self._hoppings[rn + pair] = energy

    def print_hk(self, convention: int = 1) -> None:
        """
        Print analytical Hamiltonian as the function of k-point.

        :param convention: convention for setting up the Hamiltonian
        :return: None
        """
        # Collect on-site terms
        hk = defaultdict(int)
        for i, orb in enumerate(self._orbitals):
            hk[(i, i)] = orb[1]

        # Collect hopping terms
        kpt = [sp.Symbol(_, real=True) for _ in ("ka", "kb", "kc")]
        for hop, energy in self._hoppings.items():
            rn, pair = hop[:3], hop[3:5]
            if convention == 1:
                orb_i, orb_j = pair
                pos_i = self._orbitals[orb_i][0]
                pos_j = self._orbitals[orb_j][0]
                dr = [rn[_] + pos_j[_] - pos_i[_] for _ in range(3)]
            else:
                dr = rn
            k_dot_r = kpt[0] * dr[0] + kpt[1] * dr[1] + kpt[2] * dr[2]
            phase = 2 * sp.pi * k_dot_r
            hk[pair] += energy * sp.exp(sp.I * phase)
            hk[(pair[1], pair[0])] = hk[pair].conjugate()

        # Print
        for pair, formula in hk.items():
            if pair[0] <= pair[1]:
                ham_ij = f"ham[{pair[0]}, {pair[1]}]"
                formula = sp.sympify(formula)
                print(f"{ham_ij} = {formula}")


def main():
    f = sp.Rational(1, 3)
    t = sp.Symbol("t", real=True)

    model = Model()
    model.add_orbital((f, f, 0))
    model.add_orbital((2*f, 2*f, 0))
    model.add_hopping((0, 0, 0), 0, 1, t)
    model.add_hopping((1, 0, 0), 1, 0, t)
    model.add_hopping((0, 1, 0), 1, 0, t)
    model.print_hk(1)
    model.print_hk(2)


if __name__ == "__main__":
    main()
