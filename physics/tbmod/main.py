import math

import sympy as sp
import numpy as np

from tbmod import Model


def print_hk(model: Model) -> None:
    gauges = {1: "Atomic", 2: "Lattice"}
    for conv in (1, 2):
        print(f"{gauges[conv]} gauge:")
        hk = model.get_hk(conv)
        for ii in range(model.num_orb):
            for jj in range(model.num_orb):
                print(f"ham[{ii}, {jj}] = {hk[ii, jj]}")


def kane_mele() -> None:
    # Parameters
    lat = sp.Symbol("a", real=True)
    t = sp.Symbol("t", real=True)
    lamb_so = sp.Symbol("lamb_so", real=True)
    lamb_r = sp.Symbol("lamb_r", real=True)
    lamb_nu = sp.Symbol("lamb_nu", real=True)
    f = sp.Rational(1, 3)

    # Function for generating unit vector in xOy plane
    def _ev(_theta):
        _theta *= (sp.pi / 180)
        return sp.Array([sp.cos(_theta), sp.sin(_theta), 1])

    # Lattice vectors
    vectors = sp.Matrix([_ev(-30)*lat, _ev(30)*lat, [0, 0, 1]])
    prim_cell = Model(vectors)

    # Add orbitals
    prim_cell.add_orbital((f, f, 0), energy=lamb_nu)
    prim_cell.add_orbital((f, f, 0), energy=lamb_nu)
    prim_cell.add_orbital((2*f, 2*f, 0), energy=-lamb_nu)
    prim_cell.add_orbital((2*f, 2*f, 0), energy=-lamb_nu)

    # Pauli matrices
    sigma_x = sp.Matrix([[0, 1], [1, 0]])
    sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma_z = sp.Matrix([[1, 0], [0, -1]])

    # Add 1st nearest hopping terms, which are ordinary hopping + Rashba SOC
    # dij_table is from sub-lattice A to B.
    dij_table = {(-1, 0, 0): _ev(120), (0, -1, 0): _ev(-120), (0, 0, 0): _ev(0)}
    for rn, dij in dij_table.items():
        hop_mat = t * sp.eye(2)
        hop_mat += sp.I * lamb_r * (dij[1] * sigma_x - dij[0] * sigma_y)
        for ii in range(2):
            for jj in range(2):
                prim_cell.add_hopping(rn, ii, jj+2, hop_mat[ii, jj])

    # Add 2nd nearest hopping terms, which are inter-site SOC
    # vij_table is for sub-lattice A. For B it is the opposite.
    vij_table = {(-1, 1, 0): -1, (-1, 0, 0): 1, (0, -1, 0): -1}
    for rn, vij in vij_table.items():
        hop_mat = sp.I * lamb_so * vij * sigma_z
        for ii in range(2):
            for jj in range(2):
                hij = hop_mat[ii, jj]
                prim_cell.add_hopping(rn, ii, jj, hij)
                prim_cell.add_hopping(rn, ii+2, jj+2, -hij)
    print_hk(prim_cell)


def graphene() -> None:
    # Define constants and symbols
    lattice = 2.46 * np.array([
        [1.0, 0.0, 0.0],
        [0.5, math.sqrt(3)/2, 0.0],
        [0.0, 0.0, 1.0]
    ])
    f = sp.Rational(1, 3)
    t = sp.Symbol("t", real=True)

    # Build the model
    model = Model(lattice)
    model.add_orbital((f, f, 0))
    model.add_orbital((2*f, 2*f, 0))
    model.add_hopping((0, 0, 0), 0, 1, t)
    model.add_hopping((1, 0, 0), 1, 0, t)
    model.add_hopping((0, 1, 0), 1, 0, t)

    # Print analytical Hamiltonian
    print_hk(model)

    # Plot model
    model.plot(orb_color=["r", "b"], hop_color="gray")

    # Print neighbours
    neighbors = model.find_neighbors(a_max=2, b_max=2, max_distance=2.0)
    for term in neighbors:
        print(term.rn, term.pair, term.distance)


if __name__ == "__main__":
    kane_mele()
    graphene()
