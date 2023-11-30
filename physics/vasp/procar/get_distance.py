#! /usr/bin/env python
from ase.io import read


def get_dist(atoms, index_i, index_j):
    dist = 0.0
    for j in index_j:
        dist += atoms.get_distance(index_i, j)
    dist /= len(index_j)
    return dist


def main():
    for i in ["0.98", "0.99", "1.00", "1.01", "1.02"]:
        atoms = read(f"{i}/CONTCAR.final", index=-1)

        # Cell parameters
        a, c = atoms.cell[0, 0], atoms.cell[2, 2]

        # Equatorial VO2 bond length
        v_eq = get_dist(atoms, 1, [2 ,3])

        # Axial VO2 bond length
        v_ax = get_dist(atoms, 1, [4 ,5])

        print(f"{i}, a: {a:.3f}, c: {c:.3f}, v_eq: {v_eq:.3f}, v_ax: {v_ax:.3f}")


if __name__ == "__main__":
    main()

