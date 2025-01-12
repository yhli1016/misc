import math

import numpy as np

import tbplas as tb


def add_hop(prim_cell, max_distance) -> None:
    neighbors = tb.find_neighbors(prim_cell, max_distance=max_distance)
    for term in neighbors:
        rn = term.rn
        i, j = term.pair
        prim_cell.add_hopping(rn, i, j, 1.0)


def add_hop_diag(prim_cell) -> None:
    num_orb = prim_cell.num_orb
    for i in range(num_orb//2):
        try:
            prim_cell.add_hopping((0, 0, 0), i, i+num_orb//2, 1.0)
        except:
            pass


def main():
    radius = 0.5
    num_points = 20
    lat_len = 2.0 * radius
    theta = 2 * math.pi / num_points * np.linspace(0, num_points-1, num_points)
    x = radius * np.cos(theta) + 0.5 * lat_len
    y = radius * np.sin(theta) + 0.5 * lat_len

    prim_cell = tb.PrimitiveCell(lat_vec=lat_len*np.eye(3), unit=tb.NM)
    for i in range(x.shape[0]):
        prim_cell.add_orbital_cart([x[i], y[i], 0.0], unit=tb.NM)

    add_hop(prim_cell, 0.2)
    add_hop_diag(prim_cell)
    prim_cell.plot(with_cells=False, hop_as_arrows=False)


if __name__ == "__main__":
    main()
