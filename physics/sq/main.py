import time
from itertools import product
from copy import deepcopy

import numpy as np
import scipy.linalg.lapack as lapack
import sympy as sp

from sq import SPStates, Fermion, Operator


def eval_matrix_elements(operator, basis):
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            result = operator.eval(bra, ket)
            if result != 0:
                print(f"H({ib}, {ik}) = {result}")


def get_eig_val(sp_states, num_particle, hop_terms, u_terms):
    # Set Hamiltonian matrix
    basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]
    ham = np.zeros((len(basis), len(basis)))
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            ham[ib, ik] += hop_terms.eval(bra, ket)
            ham[ib, ik] += u_terms.eval(bra, ket)

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


def split_det(det_list, i):
    result = []
    for det in det_list:
        d1, d2 = list(deepcopy(det)), list(deepcopy(det))
        if i == 0:
            d1[i] = ("u-x", d1[i][1], d1[i][2])
            d2[i] = ("u", 0, 0)
        elif i == 1:
            d1[i] = (d1[i][0], "u-x", d1[i][2])
            d2[i] = (0, "u", 0)
        else:
            d1[i] = (d1[i][0], d1[i][1], "u-x")
            d2[i] = (0, 0, "u")
        d1, d2 = tuple(d1), tuple(d2)
        result.append(d1)
        result.append(d2)
    return result


def print_det(det):
    for row in det:
        for v in row:
            print("%4s" % v, end="")
        print()
    print("------------")


def test_u():
    # Define single-particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    idx = sp_states.index

    # Define operators
    t, u = sp.symbols("t U", real=True)
    hop_terms = Operator()
    hop_terms.add_hop(idx['1+'], idx['2+'], t, True)
    hop_terms.add_hop(idx['1-'], idx['2-'], t, True)
    u_terms = Operator()
    u_terms.add_hubbard(idx['1+'], idx['1-'], u)
    u_terms.add_hubbard(idx['2+'], idx['2-'], u)

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
    # On-site terms
    hop_terms.add_ons(idx['1+'], -0.2)
    hop_terms.add_ons(idx['1-'], 0.5)
    hop_terms.add_ons(idx['2+'], -1.0)
    hop_terms.add_ons(idx['2-'], 0.1)
    hop_terms.add_ons(idx['3+'], 0.3)
    hop_terms.add_ons(idx['3-'], -0.5)
    hop_terms.add_ons(idx['4+'], 1.0)
    hop_terms.add_ons(idx['4-'], 0.1)
    # 1st nearest hopping terms
    hop_terms.add_hop(idx['1+'], idx['2+'], -1.0, True)
    hop_terms.add_hop(idx['2+'], idx['3+'], 2.5, True)
    hop_terms.add_hop(idx['3+'], idx['4+'], -1.5, True)
    hop_terms.add_hop(idx['1-'], idx['2-'], -2.0, True)
    hop_terms.add_hop(idx['2-'], idx['3-'], 1.2, True)
    hop_terms.add_hop(idx['3-'], idx['4-'], -1.0, True)
    # 2nd nearest hopping terms
    hop_terms.add_hop(idx['1+'], idx['3+'], -0.5, True)
    hop_terms.add_hop(idx['2+'], idx['4+'], 0.3, True)
    hop_terms.add_hop(idx['1-'], idx['3-'], 0.5, True)
    hop_terms.add_hop(idx['2-'], idx['4-'], -0.3, True)
    # Hubbard terms
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


def test_det():
    det = (("2u-x", "t", 0),
           ("t", "2u-x", "t"),
           (0, "t", "2u-x"))
    det_list = [det]
    det_list = split_det(det_list, 0)
    det_list = split_det(det_list, 1)
    det_list = split_det(det_list, 2)
    for det in det_list:
        print_det(det)


if __name__ == "__main__":
    test_u()
    test_eig()
    test_speed()
    test_det()