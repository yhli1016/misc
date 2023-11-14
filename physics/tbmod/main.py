import math

import sympy as sp
import numpy as np

from tbmod import Model


def main():
    lattice = 2.46 * np.array([
        [1.0, 0.0, 0.0],
        [0.5, math.sqrt(3)/2, 0.0],
        [0.0, 0.0, 1.0]
    ])
    f = sp.Rational(1, 3)
    t = sp.Symbol("t", real=True)

    model = Model(lattice)
    model.add_orbital((f, f, 0))
    model.add_orbital((2*f, 2*f, 0))
    model.add_hopping((0, 0, 0), 0, 1, t)
    model.add_hopping((1, 0, 0), 1, 0, t)
    model.add_hopping((0, 1, 0), 1, 0, t)
    model.print_hk(1)
    model.print_hk(2)

    neighbors = model.find_neighbors(a_max=2, b_max=2, max_distance=2.0)
    for term in neighbors:
        print(term.rn, term.pair, term.distance)


if __name__ == "__main__":
    main()
