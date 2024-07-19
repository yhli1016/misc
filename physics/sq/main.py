import numpy as np
import scipy.linalg.lapack as lapack
import sympy as sp
import sympy.physics.secondquant as qn

from sq import SPStates, Fermion, Operator, OperatorSympy


def eval_matrix_elements(operator, basis):
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            result = operator.eval(bra, ket)
            if result != 0:
                print(f"H({ib}, {ik}) = {result}")


def get_u_contrib(basis, hop_terms, u_terms):
    # Set Hamiltonian and U matrices
    t_mat = np.zeros((len(basis), len(basis)))
    u_mat = np.zeros((len(basis), len(basis)))
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            t_mat[ib, ik] += hop_terms.eval(bra, ket)
            u_mat[ib, ik] += u_terms.eval(bra, ket)

    # Diagonalization
    eig_val, eig_vec, info = lapack.zheev(t_mat + u_mat)
    gs = eig_vec[:, 0]
    t_contrib = np.matmul(gs.conjugate(), np.matmul(t_mat, gs))
    u_contrib = np.matmul(gs.conjugate(), np.matmul(u_mat, gs))
    print(t_contrib, u_contrib, abs(u_contrib / t_contrib))


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


def test_u_sympy():
    # Define single-particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    idx = sp_states.index_sympy

    # Define operators
    t, u = sp.symbols("t U", real=True)
    hop_terms = OperatorSympy()
    hop_terms.add_hop(idx['1+'], idx['2+'], t, True)
    hop_terms.add_hop(idx['1-'], idx['2-'], t, True)
    u_terms = OperatorSympy()
    u_terms.add_hubbard(idx['1+'], idx['1-'], u)
    u_terms.add_hubbard(idx['2+'], idx['2-'], u)

    # Define Fock states
    basis = [
            qn.FKet([idx['1+'], idx['2+']]),
            qn.FKet([idx['1-'], idx['2-']]),
            qn.FKet([idx['1+'], idx['2-']]),
            qn.FKet([idx['1-'], idx['2+']]),
            qn.FKet([idx['2+'], idx['2-']]),
            qn.FKet([idx['1+'], idx['1-']]),
        ]

    # Evaluate matrix elements
    print("Hopping terms:")
    eval_matrix_elements(hop_terms, basis)
    print("Hubbard terms:")
    eval_matrix_elements(u_terms, basis)


def test_u_contrib():
    num_site = 5
    # Define single-particle states
    sp_states = SPStates()
    for i in range(num_site):
        sp_states.append(f"{i}+")
        sp_states.append(f"{i}-")
    idx = sp_states.index

    # Define operators
    t = 1.0
    u = 2.0
    hop_terms = Operator()
    for i in range(num_site - 1):
        hop_terms.add_hop(idx[f"{i}+"], idx[f"{i+1}+"], t, True)
        hop_terms.add_hop(idx[f"{i}-"], idx[f"{i+1}-"], t, True)
    u_terms = Operator()
    for i in range(num_site):
        u_terms.add_hubbard(idx[f"{i}+"], idx[f"{i}-"], u)

    # Define Fock states
    for num_particle in range(1, num_site * 2+1):
        basis = [Fermion(_) for _ in sp_states.combinations(num_particle)]
        get_u_contrib(basis, hop_terms, u_terms)


if __name__ == "__main__":
    test_u()
    test_u_sympy()
    test_u_contrib()
