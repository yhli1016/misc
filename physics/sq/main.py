import numpy as np
import scipy.linalg.lapack as lapack

from sq import ReferenceStates, FermionState, TwoBody, Hubbard


def test_u():
    # Define reference states
    ref_states = ReferenceStates()
    ref_states.add_state('+1')
    ref_states.add_state('-1')
    ref_states.add_state('+2')
    ref_states.add_state('-2')
    num_ref = ref_states.size

    # Define Fock states
    def f(i):
        return ref_states.index(i)
    basis = [
            FermionState(num_ref, [f('+1'), f('+2')]),
            FermionState(num_ref, [f('-1'), f('-2')]),
            FermionState(num_ref, [f('+1'), f('-2')]),
            FermionState(num_ref, [f('-1'), f('+2')]),
            FermionState(num_ref, [f('+2'), f('-2')]),
            FermionState(num_ref, [f('+1'), f('-1')]),
        ]

    # Define operators
    hop_terms = [
        TwoBody(f('+1'), f('+2')),
        TwoBody(f('-1'), f('-2')),
        TwoBody(f('+2'), f('+1')),
        TwoBody(f('-2'), f('-1')),
    ]
    u_terms = [
        Hubbard(f('+1'), f('-1')),
        Hubbard(f('+2'), f('-2')),
    ]

    # Evaluate matrix elements
    print("Hopping terms:")
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in hop_terms:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    print(ib, ik, operator.indices, prod)

    print("Hubbard terms:")
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in u_terms:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    print(ib, ik, operator.indices, prod)


def test_eig():
    # Define reference states
    ref_states = ReferenceStates()
    ref_states.add_state('+1')
    ref_states.add_state('-1')
    ref_states.add_state('+2')
    ref_states.add_state('-2')
    ref_states.add_state('+3')
    ref_states.add_state('-3')
    ref_states.add_state('+4')
    ref_states.add_state('-4')
    num_ref = ref_states.size

    # Define Fock states
    def f(i):
        return ref_states.index(i)
    num_particle = 2
    basis = [FermionState(num_ref, _) for _ in ref_states.combinations(num_particle)]

    # Define operators
    hop_terms = [
        TwoBody(f('+1'), f('+2')),
        TwoBody(f('+2'), f('+1')),
        TwoBody(f('+2'), f('+3')),
        TwoBody(f('+3'), f('+2')),
        TwoBody(f('+3'), f('+4')),
        TwoBody(f('+4'), f('+3')),
        TwoBody(f('-1'), f('-2')),
        TwoBody(f('-2'), f('-1')),
        TwoBody(f('-2'), f('-3')),
        TwoBody(f('-3'), f('-2')),
        TwoBody(f('-3'), f('-4')),
        TwoBody(f('-4'), f('-3')),
    ]
    u_terms = [
        Hubbard(f('+1'), f('-1')),
        Hubbard(f('+2'), f('-2')),
        Hubbard(f('+3'), f('-3')),
        Hubbard(f('+4'), f('-4')),
    ]

    # Get eigenvalues and eigenstates
    basis_size = len(basis)
    a = np.zeros((basis_size, basis_size))
    t = 1.0
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in hop_terms:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    a[ib, ik] = t * prod
    eig_val, eig_vec, info = lapack.zheev(a)
    print(eig_val)


if __name__ == "__main__":
    test_eig()
