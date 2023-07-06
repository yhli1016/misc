from copy import deepcopy
from itertools import permutations, combinations
from collections import defaultdict
from typing import List, Dict, Hashable, Iterable


class SPStates:
    """
    Container for holding single-particle states.

    Attributes
    ----------
    _states: List[Hashable]
        identifiers of single particle states
    """
    def __init__(self) -> None:
        self._states = []

    def append(self, state: Hashable) -> None:
        """
        Append a new state to the list of states.

        :param state: identifier of the new state
        :return: None
        """
        if state not in self._states:
            self._states.append(state)

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

        :return: indexing array
        """
        return {state: idx for idx, state in enumerate(self._states)}


class Boson:
    """
    Fock state for Bosons.

    Attributes
    ----------
    _is_null: boolean
        whether the state is a null state
    _sign: int
        sign of the state, effective when evaluating inner products
        For the null state, we set it to 0 for safety.
    _occ: Dict[int, int]
        occupation numbers on the single particle states
    """
    def __init__(self, occupied_states: Iterable[int] = None) -> None:
        self._is_null = False
        self._sign = 1
        self._occ = defaultdict(int)
        for idx in occupied_states:
            self._occ[idx] += 1

    def __eq__(self, other) -> bool:
        """
        Check if this state equals to the other state.

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
        self._is_null = True
        self._sign = 0
        self._occ = defaultdict(int)

    def purge(self) -> None:
        """
        Remove states with zero occupation.

        :return: None
        """
        keys = list(self._occ.keys())
        for idx in keys:
            if self._occ[idx] <= 0:
                self._occ.pop(idx)

    def inner_prod(self, ket) -> int:
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
                prod = self.sign * ket.sign
            else:
                prod = 0
        return prod

    @property
    def is_null(self) -> bool:
        """Interface for the '_is_null' attribute."""
        return self._is_null

    @property
    def sign(self) -> int:
        """Interface for the '_sign' attribute."""
        return self._sign

    @property
    def occ(self) -> Dict[int, int]:
        """Interface for the '_occ' attribute."""
        return self._occ

    @property
    def is_empty(self) -> bool:
        """Check if the occupation numbers are empty."""
        return sum([abs(_) for _ in self._occ.values()]) == 0

    @property
    def is_vac(self) -> bool:
        """Check if the state is the vacuum state."""
        return not self.is_null and self.is_empty

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
                self._sign *= self.sign_factor(idx)

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
                self._sign *= self.sign_factor(idx)

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


class TwoBody:
    """Operator of c_{i+} c_j."""
    def __init__(self, i: int, j: int) -> None:
        self.indices = (i, j)

    def eval(self, bra: Boson, ket: Boson) -> int:
        ket_c = deepcopy(ket)
        ket_c.destroy(self.indices[1])
        ket_c.create(self.indices[0])
        return bra.inner_prod(ket_c)


class OnSite(TwoBody):
    """Operator of c_{i+} c_i."""
    def __init__(self, i: int) -> None:
        super().__init__(i, i)


class FourBody:
    """Operator of c_{i+} c_{j+} c_n c_m."""
    def __init__(self, i: int, j: int, n: int, m: int):
        self.indices = (i, j, n, m)

    def eval(self, bra: Boson, ket: Boson) -> int:
        ket_c = deepcopy(ket)
        ket_c.destroy(self.indices[3])
        ket_c.destroy(self.indices[2])
        ket_c.create(self.indices[1])
        ket_c.create(self.indices[0])
        return bra.inner_prod(ket_c)


class Hubbard(FourBody):
    """Operator of c_{i+} c_i c_{j+} c_j."""
    def __init__(self, i: int, j: int):
        super().__init__(i, j, j, i)
