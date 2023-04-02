from copy import deepcopy
from itertools import permutations, combinations
from abc import ABC, abstractmethod
from typing import List, Any, Iterable


class RefStates:
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
    is_null: boolean
        whether the state is a null state
    sign: int
        sign of the state, effective when evaluating inner products
        For Bosons it should always be 1, while for Fermions it can be
        either -1 or 1. For the null state, we set it to 0 for safety.
    occupations: Dict[int, int]
        occupation numbers on the single particle states
    """
    def __init__(self, occupied_states: Iterable[int] = None) -> None:
        self.is_null = False
        self.sign = 1
        self.occupations = dict()
        for idx in occupied_states:
            self.occupations[idx] = self.get_occ(idx) + 1

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
            self.purge()
            other.purge()
            is_equal = self.occupations == other.occupations
        return is_equal

    def get_occ(self, idx: int) -> int:
        """Safe wrapper to extract occupation number."""
        try:
            occ = self.occupations[idx]
        except KeyError:
            occ = 0
        return occ

    def nullify(self) -> None:
        """Set the state to the null state."""
        self.is_null = True
        self.sign = 0
        self.occupations = dict()

    def purge(self) -> None:
        """Remove states with zero occupation."""
        new_occ = {idx2: occ
                   for idx2, occ in self.occupations.items()
                   if occ > 0}
        self.occupations = new_occ

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


class Boson(FockState):
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
            occ = self.get_occ(idx)
            if occ <= 0:
                self.nullify()
            else:
                self.occupations[idx] = occ - 1

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
            self.occupations[idx] = self.get_occ(idx) + 1


class Fermion(FockState):
    def __init__(self, occupied_states: Iterable[int] = None) -> None:
        super().__init__(set(occupied_states))

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
            occ = self.get_occ(idx)
            if occ <= 0:
                self.nullify()
            else:
                self.occupations[idx] = occ - 1
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
            occ = self.get_occ(idx)
            if occ >= 1:
                self.nullify()
            else:
                self.occupations[idx] = occ + 1
                self.sign *= self.sign_factor(idx)

    def sign_factor(self, idx: int) -> int:
        """
        Get the sign factor when creating and destroying particles.

        :param int idx: index of the state on which the particle is created
            or destroyed
        :return: sign factor
        """
        power = 0
        for idx2, occ in self.occupations.items():
            if idx2 < idx and occ > 0:
                power += 1
        factor = (-1)**power
        return factor


class TwoBody:
    """Operator of c_{i+} c_j."""
    def __init__(self, i: int, j: int) -> None:
        self.indices = (i, j)

    def eval(self, bra: FockState, ket: FockState) -> int:
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

    def eval(self, bra: FockState, ket: FockState) -> int:
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
