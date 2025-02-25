from copy import deepcopy
from itertools import permutations, combinations
from collections import defaultdict
from typing import Dict, Hashable, Iterable, Tuple, Union

import sympy as sp
import sympy.physics.secondquant as qn


c_type = Union[int, float, complex, sp.Basic]
term_type = Union[Tuple[int, int], Tuple[int, int, int, int]]


class SPStates:
    """
    Container for holding single-particle states.

    Attributes
    ----------
    _states: Dict[Hashable, int]
        tags and indices of single particle states, with indices starting from 0
    """
    def __init__(self) -> None:
        self._states = dict()

    def append(self, tag: Hashable) -> None:
        """
        Append a new state to the list of states.

        :param tag: tag of the new state
        :return: None
        :raises RuntimeError: if the new state is duplicate
        """
        if tag not in self._states.keys():
            num_states = len(self._states)
            self._states[tag] = num_states
        else:
            raise RuntimeError(f"Duplicate state {tag}")

    def permutations(self, num_states: int = 1) -> permutations:
        """
        Get the permutations of state indices, for constructing operators.

        :param num_states: number of selected states
        :return: all possible permutations
        """
        return permutations(range(self.size), num_states)

    def combinations(self, num_states: int = 1) -> combinations:
        """
        Get the combinations of state indices, for constructing Fock states.

        :param num_states: number of selected states
        :return: all possible combinations
        """
        return combinations(range(self.size), num_states)

    def permutations_sympy(self, num_states: int = 1) -> permutations:
        """
        Get the permutations of state indices, for constructing operators.

        The indices start from 1, which is suitable for sympy.

        :param num_states: number of selected states
        :return: all possible permutations
        """
        return permutations(range(1, self.size+1), num_states)

    def combinations_sympy(self, num_states: int = 1) -> combinations:
        """
        Get the combinations of state indices, for constructing Fock states.

        The indices start from 1, which is suitable for sympy.

        :param num_states: number of selected states
        :return: all possible combinations
        """
        return combinations(range(1, self.size+1), num_states)

    @property
    def size(self) -> int:
        """
        Get the number of single particle states.

        :return: number of single particle states
        """
        return len(self._states)

    @property
    def index(self) -> Dict[Hashable, int]:
        """
        Get the dictionary for indexing single particle states.

        :return: indexing dictionary
        """
        return self._states

    @property
    def index_sympy(self) -> Dict[Hashable, int]:
        """
        Get the dictionary for indexing single particle states, but starting
        from 1, which is required by the second quantization module of sympy.

        :return: indexing dictionary
        """
        idx = dict()
        for key, value in self._states.items():
            idx[key] = value + 1
        return idx


class Boson:
    """
    Fock state for Bosons.

    Attributes
    ----------
    _coeff: c_type
        coefficient of the state, effective when evaluating inner products
        For the null state, we set it to 0 for safety.
    _occ: Dict[int, int] or None
        occupation numbers on the single particle states
        None for a null state (numerically zero in second quantization regime).
    """
    def __init__(self, occupied_states: Iterable[int] = None) -> None:
        self._coeff = 1
        self._occ = defaultdict(int)
        for idx in occupied_states:
            self._occ[idx] += 1

    def __eq__(self, other) -> bool:
        """
        Check if this state equals to the other state by comparing their
        occupation numbers.

        :param other: other state to compare
        :return: whether this state equals to the other state
        """
        # Null state equals to itself only. For vacuum and common states,
        # equality is determined from the occupation numbers.
        if self.is_null and other.is_null:
            is_equal = True
        elif self.is_null and not other.is_null:
            is_equal = False
        elif not self.is_null and other.is_null:
            is_equal = False
        else:
            self.purge()
            other.purge()
            is_equal = self.occ == other.occ
        return is_equal

    def nullify(self) -> None:
        """
        Set the state to the null state.

        :return: None
        """
        self._coeff = 0
        self._occ = None

    def purge(self) -> None:
        """
        Remove states with zero occupation.

        :return: None
        """
        keys = list(self._occ.keys())
        for idx in keys:
            if self._occ[idx] <= 0:
                self._occ.pop(idx)

    def inner_prod(self, ket) -> c_type:
        """
        Evaluate the inner product <self|ket>.

        :param ket: the right operand
        :return: the inner product
        """
        # Inner product involving null state is always zero.
        # Otherwise, check if the states are equal.
        if self.is_null or ket.is_null:
            prod = 0
        else:
            if self == ket:
                prod = self.coeff * ket.coeff
            else:
                prod = 0
        return prod

    @property
    def coeff(self) -> c_type:
        """Interface for the '_coeff' attribute."""
        return self._coeff

    @property
    def occ(self) -> Dict[int, int]:
        """Interface for the '_occ' attribute."""
        return self._occ

    @property
    def is_null(self) -> bool:
        """Check if the state is the null state."""
        return self._occ is None

    @property
    def is_vac(self) -> bool:
        """Check if the state is the vacuum state."""
        sum_occ = sum([abs(_) for _ in self._occ.values()])
        return (self._occ is not None) and (sum_occ == 0)

    def destroy(self, idx: int) -> None:
        """
        Destroy a particle on given single particle state.

        :param int idx: index of the state
        :return: None.
        """
        # Destroying a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            occ = self._occ[idx]
            if occ <= 0:
                self.nullify()
            else:
                self._coeff *= sp.sqrt(occ)
                self._occ[idx] = occ - 1

    def create(self, idx: int) -> None:
        """
        Create a particle given single particle state.

        :param int idx: index of the state
        :return: None.
        """
        # Creating a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            occ = self._occ[idx]
            self._coeff *= sp.sqrt(occ + 1)
            self._occ[idx] += 1


class Fermion(Boson):
    def __init__(self, occupied_states: Iterable[int] = None) -> None:
        super().__init__(set(occupied_states))

    def destroy(self, idx: int) -> None:
        """
        Destroy a particle on given single particle state.

        :param int idx: index of the state
        :return: None.
        """
        # Destroying a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            occ = self._occ[idx]
            if occ <= 0:
                self.nullify()
            else:
                self._occ[idx] = occ - 1
                self._coeff *= self.sign_factor(idx)

    def create(self, idx: int) -> None:
        """
        Create a particle given single particle state.

        :param int idx: index of the state
        :return: None.
        """
        # Creating a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            occ = self._occ[idx]
            if occ >= 1:
                self.nullify()
            else:
                self._occ[idx] = occ + 1
                self._coeff *= self.sign_factor(idx)

    def sign_factor(self, idx: int) -> int:
        """
        Get the sign factor when creating and destroying particles.

        :param int idx: index of the state
        :return: sign factor
        """
        power = 0
        for idx2, occ in self._occ.items():
            if idx2 < idx and occ > 0:
                power += 1
        factor = (-1)**power
        return factor


class Operator:
    """
    Class for representing an operator consisted of two-body and four-body terms
    in second quantized form.

    Attributes
    ----------
    _terms: Dict[term_type, c_type]
        indices and coefficients of the terms, with keys being tuples of 2 or 4
        integers
    """
    def __init__(self) -> None:
        self._terms = dict()

    def __add_term(self, term: term_type, coeff: c_type = 0) -> None:
        """Add a new term to the operator."""
        if term not in self._terms.keys():
            self._terms[term] = coeff
        else:
            raise ValueError(f"Duplicate term {term}")

    def add_2bd(self, i: int, j: int, coeff: c_type = 0) -> None:
        """Add a general two-body term c_{i+} c_j."""
        self.__add_term((i, j), coeff)

    def add_4bd(self, i: int, j: int, n: int, m: int,
                coeff: c_type = 0) -> None:
        """Add a general four-body term c_{i+} c_{j+} c_n c_m."""
        self.__add_term((i, j, n, m), coeff)

    def add_ons(self, i: int, coeff: c_type = 0) -> None:
        """Add an on-site term c_{i+} c_i."""
        self.add_2bd(i, i, coeff)

    def add_hop(self, i: int, j: int, coeff: c_type = 0,
                with_conj: bool = False) -> None:
        """Add a hopping term c_{i+} c_j with i != j."""
        if i != j:
            self.add_2bd(i, j, coeff)
            if with_conj:
                self.add_2bd(j, i, coeff.conjugate())
        else:
            raise ValueError(f"Hopping term require {i} != {j}")

    def add_hubbard(self, i: int, j: int, coeff: c_type = 0) -> None:
        """
        Add a Hubbard term c_{i+} c_i c_{j+} c_j.

        The derivation follows
            (i+, i, j+, j) = -(i+, j+, i, j) = (i+, j+, j, i)
        which requires i != j.
        """
        if i != j:
            self.add_4bd(i, j, j, i, coeff)
        else:
            raise ValueError(f"Hubbard term require {i} != {j}")

    def eval(self, bra: Boson, ket: Boson) -> c_type:
        """
        Evaluate the matrix element between two Fock states.

        :param bra: left operand
        :param ket: right operand
        :return: the matrix element
        """
        result = 0
        for term, coeff in self._terms.items():
            ket_c = deepcopy(ket)
            if len(term) == 2:
                ket_c.destroy(term[1])
                ket_c.create(term[0])
            else:
                ket_c.destroy(term[3])
                ket_c.destroy(term[2])
                ket_c.create(term[1])
                ket_c.create(term[0])
            prod = bra.inner_prod(ket_c)
            result += coeff * prod
        return result


class OperatorSympy(Operator):
    """
    Class for representing an operator consisted of two-body and four-body terms
    in second quantized form for Fermions.

    Attributes
    ----------
    _terms: sp.Add
        terms in the operator
    """
    def __init__(self) -> None:
        super().__init__()
        self._terms = 0

    def add_2bd(self, i: int, j: int, coeff: c_type = 0) -> None:
        """Add a general two-body term c_{i+} c_j."""
        self._terms += coeff * qn.Fd(i) * qn.F(j)

    def add_4bd(self, i: int, j: int, n: int, m: int,
                coeff: c_type = 0) -> None:
        """Add a general four-body term c_{i+} c_{j+} c_n c_m."""
        self._terms += coeff * qn.Fd(i) * qn.Fd(j) * qn.F(n) * qn.F(m)

    def eval(self, bra: qn.FermionState, ket: qn.FermionState) -> c_type:
        """
        Evaluate the matrix element between two Fock states.

        :param bra: left operand
        :param ket: right operand
        :return: the matrix element
        """
        result = qn.apply_operators(qn.Dagger(bra) * self._terms * ket)
        return result
