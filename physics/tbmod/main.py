import math

import sympy as sp
import numpy as np

from tbmod import Model


def main():
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
    hk1 = model.get_hk(1)
    hk2 = model.get_hk(2)
    for i in range(model.num_orb):
        for j in range(model.num_orb):
            print(f"h[{i}, {j}] = {hk1[i, j]}")
    for i in range(model.num_orb):
        for j in range(model.num_orb):
            print(f"h[{i}, {j}] = {hk2[i, j]}")

    # Plot model
    model.plot(orb_color=["r", "b"], hop_color="gray")

    # Print neighbours
    neighbors = model.find_neighbors(a_max=2, b_max=2, max_distance=2.0)
    for term in neighbors:
        print(term.rn, term.pair, term.distance)


if __name__ == "__main__":
    main()
