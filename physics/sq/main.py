import time
from itertools import product

import numpy as np
import scipy.linalg.lapack as lapack

from sq import SPStates, Fermion, Operator


def eval_matrix_elements(operator, basis):
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            result = operator.eval(bra, ket)
            if len(result) != 0:
                print(f"H({ib}, {ik}) =")
                for term, coeff in result.items():
                    print(f"\t{coeff} * h{term}")


def get_eig_val(sp_states, num_particle, hop_terms, u_terms):
    # Set Hamiltonian matrix
    basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]
    ham = np.zeros((len(basis), len(basis)))
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            # Hopping terms
            result = hop_terms.eval(bra, ket)
            ham[ib, ik] += sum(result.values())

            # Hubbard term
            result = u_terms.eval(bra, ket)
            ham[ib, ik] += sum(result.values())

    # Diagonalization
    eig_val, eig_vec, info = lapack.zheev(ham)
    return eig_val


def add(a1, a2):
    result = []
    for i in product(a1, a2):
        result.append(i[0] + i[1])
    return result


def diff(a1, a2):
    for x in a1:
        x_in_a2 = False
        for y in a2:
            if abs(y - x) < 1.0e-7:
                x_in_a2 = True
                break
        if not x_in_a2:
            print(x)


def test_u():
    # Define single-particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    idx = sp_states.index

    # Define operators
    hop_terms = Operator()
    hop_terms.add_hop(idx['1+'], idx['2+'], 1, True)
    hop_terms.add_hop(idx['1-'], idx['2-'], 1, True)
    u_terms = Operator()
    u_terms.add_hubbard(idx['1+'], idx['1-'], 1)
    u_terms.add_hubbard(idx['2+'], idx['2-'], 1)

    # Define Fock states
    basis = [
            Fermion([idx['1+'], idx['2+']]),
            Fermion([idx['1-'], idx['2-']]),
            Fermion([idx['1+'], idx['2-']]),
            Fermion([idx['1-'], idx['2+']]),
            Fermion([idx['2+'], idx['2-']]),
            Fermion([idx['1+'], idx['1-']]),
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

    # Define operators
    hop_terms = Operator()
    hop_terms.add_ons(idx['1+'], -0.2)
    hop_terms.add_ons(idx['1-'], 0.5)
    hop_terms.add_ons(idx['2+'], -1.0)
    hop_terms.add_ons(idx['2-'], 0.1)
    hop_terms.add_ons(idx['3+'], 0.3)
    hop_terms.add_ons(idx['3-'], -0.5)
    hop_terms.add_ons(idx['4+'], 1.0)
    hop_terms.add_ons(idx['4-'], 0.1)
    hop_terms.add_hop(idx['1+'], idx['2+'], -1.0, True)
    hop_terms.add_hop(idx['2+'], idx['3+'], 2.5, True)
    hop_terms.add_hop(idx['3+'], idx['4+'], -1.5, True)
    hop_terms.add_hop(idx['1-'], idx['2-'], -2.0, True)
    hop_terms.add_hop(idx['2-'], idx['3-'], 1.2, True)
    hop_terms.add_hop(idx['3-'], idx['4-'], -1.0, True)
    u_terms = Operator()
    u_terms.add_hubbard(idx['1+'], idx['1-'], 0.0)
    u_terms.add_hubbard(idx['2+'], idx['2-'], 0.0)
    u_terms.add_hubbard(idx['3+'], idx['3-'], 0.0)
    u_terms.add_hubbard(idx['4+'], idx['4-'], 0.0)

    # Actual eigenvalues
    eig1 = get_eig_val(sp_states, 1, hop_terms, u_terms)
    eig2 = get_eig_val(sp_states, 2, hop_terms, u_terms)
    eig3 = get_eig_val(sp_states, 3, hop_terms, u_terms)
    eig4 = get_eig_val(sp_states, 4, hop_terms, u_terms)

    # Sum of eigenvalues
    eig11 = add(eig1, eig1)
    eig12 = add(eig1, eig2)
    eig13 = add(eig1, eig3)
    eig22 = add(eig2, eig2)

    # Check
    print("\nTwo particle")
    diff(eig2, eig11)
    print("\nTree particle")
    diff(eig3, eig12)
    print("\nFour particle 13")
    diff(eig4, eig13)
    print("\nFour particle 22")
    diff(eig4, eig22)


def test_speed():
    # Define single particle states
    sp_states = SPStates()
    for i in range(10):
        sp_states.append(i)

    # Define Fock states
    num_particle = 2
    basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]

    # Define operators
    hop_terms = Operator()
    for i, j in sp_states.permutations(2):
        hop_terms.add_hop(i, j)

    # Benchmark
    t0 = time.time()
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            hop_terms.eval(bra, ket)
    t1 = time.time()
    print(t1 - t0)


if __name__ == "__main__":
    test_u()
    test_eig()
    test_speed()
