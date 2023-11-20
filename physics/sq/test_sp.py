import sympy as sp
import sympy.physics.secondquant as qn

from sq import SPStates


def hop(i: int, j: int) -> sp.Mul:
    return qn.Fd(i) * qn.F(j)


def hubbard(i: int, j: int) -> sp.Mul:
    return qn.Fd(i) * qn.F(i) * qn.Fd(j) * qn.F(j)


def test_u2():
    # Define single-particle states
    sp_states = SPStates()
    sp_states.append('1+')
    sp_states.append('1-')
    sp_states.append('2+')
    sp_states.append('2-')
    idx = sp_states.index_sympy

    # Define operators
    t, u = sp.symbols("t U", real=True)
    hop_terms = t * (hop(idx['1+'], idx['2+']) + hop(idx['1-'], idx['2-']))
    hop_terms += t * (hop(idx['2+'], idx['1+']) + hop(idx['2-'], idx['1-']))
    u_terms = u * (hubbard(idx['1+'], idx['1-']) + hubbard(idx['2+'],idx['2-']))
    hop_terms += u_terms

    # Define Fock states
    basis = [
            qn.FKet([idx['1+'], idx['2+']]),
            qn.FKet([idx['1-'], idx['2-']]),
            qn.FKet([idx['1+'], idx['2-']]),
            qn.FKet([idx['1-'], idx['2+']]),
            qn.FKet([idx['2+'], idx['2-']]),
            qn.FKet([idx['1+'], idx['1-']]),
        ]

    ham = sp.zeros(len(basis))
    for ib, bra in enumerate(basis):
        for ik, ket in enumerate(basis):
            ham[ib, ik] = qn.apply_operators(qn.Dagger(bra) * hop_terms * ket)
    print(ham)


if __name__ == "__main__":
    test_u2()
