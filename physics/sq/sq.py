from copy import deepcopy
from itertools import permutations, combinations
from abc import ABC, abstractmethod
from typing import List, Any


class ReferenceStates:
    """
    Container for holding single-particle reference states.

    Attributes
    ----------
    states: List[Any]
        identifiers of single particle states
    """
    def __init__(self) -> None:
        self.states = []

    def add_state(self, state: Any) -> None:
        if state not in self.states:
            self.states.append(state)

    def get_state(self, idx: int) -> Any:
        return self.states[idx]

    def index(self, state: Any) -> int:
        return self.states.index(state)

    def permutations(self, num_states: int = 1) -> permutations:
        return permutations(range(self.size), num_states)

    def combinations(self, num_states: int = 1) -> combinations:
        return combinations(range(self.size), num_states)

    @property
    def size(self) -> int:
        return len(self.states)


class FockState(ABC):
    """
    Base class of Fock states for Bosons and Fermions.

    Attributes
    ----------
    num_ref: int
        number of reference single particle states
    is_null: boolean
        whether the state is a null state
    sign: int
        sign of the state, effective when evaluating inner products
    occupations: List[int]
        occupation numbers on the single particle states
    """
    def __init__(self, num_ref: int = 1,
                 occupied_states: List[int] = None) -> None:
        self.num_ref = num_ref
        self.is_null = False
        self.sign = 1
        self.occupations = [0 for _ in range(self.num_ref)]
        for idx in occupied_states:
            self.occupations[idx] += 1

    def __eq__(self, other) -> bool:
        """
        Check if this state equals to the other state.

        :param FermiState other: other state to compare
        :return: whether this state equals to the other state
        """
        # Null state equals to itself only.
        if self.is_null and other.is_null:
            is_equal = True
        elif self.is_null and not other.is_null:
            is_equal = False
        elif not self.is_null and other.is_null:
            is_equal = False
        #  For other states (vacuum and common), equality is determined
        #  from the occupation numbers.
        else:
            is_equal = self.occupations == other.occupations
        return is_equal

    def nullify(self) -> None:
        """
        Set the state to null state.

        :return: None.
        """
        self.is_null = True
        self.sign = 0
        self.occupations = [0 for _ in range(self.num_ref)]

    @abstractmethod
    def destroy(self, idx: int) -> None:
        pass

    @abstractmethod
    def create(self, idx: int) -> None:
        pass

    def inner_prod(self, ket) -> int:
        """
        Evaluate the inner product <self|ket>.

        :param FermiState ket: the right operand
        :return: the inner product
        """
        # Inner product involving null state is always zero.
        if self.is_null or ket.is_null:
            prod = 0
        # Otherwise, check if the states are equal.
        else:
            if self == ket:
                prod = self.sign * ket.sign
            else:
                prod = 0
        return prod

    @property
    def is_empty(self) -> bool:
        """Check if the occupation numbers are empty."""
        return sum([abs(_) for _ in self.occupations]) == 0

    @property
    def is_vac(self) -> bool:
        """Check if the state is the vacuum state."""
        return not self.is_null and self.is_empty


class BosonState(FockState):
    def destroy(self, idx: int) -> None:
        """
        Destroy a particle on given single particle state.

        :param int idx: index of the state on which the particle is destroyed
        :return: None.
        """
        # Destroying a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            if self.occupations[idx] <= 0:
                self.nullify()
            else:
                self.occupations[idx] -= 1

    def create(self, idx: int) -> None:
        """
        Create a particle given single particle state.

        :param int idx: index of the state on which the particle is created
        :return: None.
        """
        # Creating a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            self.occupations[idx] += 1


class FermionState(FockState):
    def __init__(self, num_ref: int,
                 occupied_states: List[int] = None) -> None:
        super().__init__(num_ref, occupied_states)
        for i in self.occupations:
            if i not in (0, 1):
                raise ValueError(f"Occupation numbers of Fermions should be"
                                 f" either 0 or 1")

    def destroy(self, idx: int) -> None:
        """
        Destroy a particle on given single particle state.

        :param int idx: index of the state on which the particle is destroyed
        :return: None.
        """
        # Destroying a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            if self.occupations[idx] <= 0:
                self.nullify()
            else:
                self.occupations[idx] -= 1
                self.sign *= self.sign_factor(idx)

    def create(self, idx: int) -> None:
        """
        Create a particle given single particle state.

        :param int idx: index of the state on which the particle is created
        :return: None.
        """
        # Creating a particle on the null state always yields itself.
        if self.is_null:
            pass
        else:
            if self.occupations[idx] >= 1:
                self.nullify()
            else:
                self.occupations[idx] += 1
                self.sign *= self.sign_factor(idx)

    def sign_factor(self, idx: int) -> int:
        """
        Get the sign factor when creating and destroying particles.

        :param int idx: index of the state on which the particle is created
            or destroyed
        :return: sign factor
        """
        power = 0
        for j in self.occupations[:idx]:
            if j > 0:
                power += 1
        factor = (-1)**power
        return factor


class TwoBody:
    """Matrix element of c_{i+} c_j."""
    def __init__(self, i: int, j: int) -> None:
        self.indices = (i, j)

    def eval(self, bra: FockState, ket: FockState) -> int:
        bra_c = deepcopy(bra)
        ket_c = deepcopy(ket)
        bra_c.destroy(self.indices[0])
        ket_c.destroy(self.indices[1])
        return bra_c.inner_prod(ket_c)


class OnSite(TwoBody):
    """Matrix element of c_{i+} c_i."""
    def __init__(self, i: int) -> None:
        super().__init__(i, i)


class FourBody:
    """Matrix element of c_{i+} c_j c_{k+} c_l."""
    def __init__(self, i: int, j: int, m: int, n: int):
        self.indices = (i, j, m, n)

    def eval(self, bra: FockState, ket: FockState) -> int:
        bra_c = deepcopy(bra)
        ket_c = deepcopy(ket)
        bra_c.destroy(self.indices[0])
        bra_c.create(self.indices[1])
        ket_c.destroy(self.indices[3])
        ket_c.create(self.indices[2])
        return bra_c.inner_prod(ket_c)


class Hubbard(FourBody):
    """Matrix element of c_{i+} c_i c_{j+} c_j."""
    def __init__(self, i: int, j: int):
        super().__init__(i, i, j, j)
