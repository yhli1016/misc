import time

import numpy as np
import scipy.linalg.lapack as lapack

from sq import SPStates, Fermion, TwoBody, Hubbard


def eval_matrix_elements(operators, basis):
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in operators:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    print(ib, ik, operator.indices, prod)


def test_u():
    # Define single-particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    idx = sp_states.index

    # Define Fock states
    basis = [
            Fermion([idx['1+'], idx['2+']]),
            Fermion([idx['1-'], idx['2-']]),
            Fermion([idx['1+'], idx['2-']]),
            Fermion([idx['1-'], idx['2+']]),
            Fermion([idx['2+'], idx['2-']]),
            Fermion([idx['1+'], idx['1-']]),
        ]

    # Define operators
    hop_terms = [
        TwoBody(idx['1+'], idx['2+']),
        TwoBody(idx['1-'], idx['2-']),
        TwoBody(idx['2+'], idx['1+']),
        TwoBody(idx['2-'], idx['1-']),
    ]
    u_terms = [
        Hubbard(idx['1+'], idx['1-']),
        Hubbard(idx['2+'], idx['2-']),
    ]

    # Evaluate matrix elements
    print("Hopping terms:")
    eval_matrix_elements(hop_terms, basis)

    print("Hubbard terms:")
    eval_matrix_elements(u_terms, basis)


def test_eig():
    # Define single particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    sp_states.append('3+')
    sp_states.append('3-')
    sp_states.append('4+')
    sp_states.append('4-')
    idx = sp_states.index

    # Define Fock states
    num_particle = 2
    basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]

    # Define operators
    hop_terms = [
        TwoBody(idx['1+'], idx['2+']),
        TwoBody(idx['2+'], idx['1+']),
        TwoBody(idx['2+'], idx['3+']),
        TwoBody(idx['3+'], idx['2+']),
        TwoBody(idx['3+'], idx['4+']),
        TwoBody(idx['4+'], idx['3+']),
        TwoBody(idx['1-'], idx['2-']),
        TwoBody(idx['2-'], idx['1-']),
        TwoBody(idx['2-'], idx['3-']),
        TwoBody(idx['3-'], idx['2-']),
        TwoBody(idx['3-'], idx['4-']),
        TwoBody(idx['4-'], idx['3-']),
    ]
    u_terms = [
        Hubbard(idx['1+'], idx['1-']),
        Hubbard(idx['2+'], idx['2-']),
        Hubbard(idx['3+'], idx['3-']),
        Hubbard(idx['4+'], idx['4-']),
    ]

    # Get eigenvalues and eigenstates
    basis_size = len(basis)
    a = np.zeros((basis_size, basis_size))
    t = 1.0
    u = 0.5
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in hop_terms:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    a[ib, ik] += t * prod
            for operator in u_terms:
                prod = operator.eval(bra, ket)
                if prod != 0:
                    a[ib, ik] += u * prod
    eig_val, eig_vec, info = lapack.zheev(a)
    print(eig_val)


def test_speed():
    # Define single particle states
    sp_states = SPStates()
    for i in range(10):
        sp_states.append(i)

    # Define Fock states
    num_particle = 2
    basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]

    # Define operators
    hop_terms = [TwoBody(i, j) for i, j in sp_states.permutations(2)]

    # Benchmark
    t0 = time.time()
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            for operator in hop_terms:
                operator.eval(bra, ket)
    t1 = time.time()
    print(t1 - t0)


if __name__ == "__main__":
    test_u()
    test_eig()
    test_speed()
