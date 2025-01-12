#! /usr/bin/env python

import numpy as np

import tbplas as tb


def main():
    lat_vec = np.eye(3, dtype=np.float64)
    prim_cell = tb.PrimitiveCell(lat_vec)
    prim_cell.add_orbital((0.25, 0.25), 1.0)
    prim_cell.add_orbital((0.75, 0.75), 2.0)
    prim_cell.add_orbital((0.25, 0.75), 3.0)
    prim_cell.add_orbital((0.75, 0.25), 4.0)
    for orb_i in range(prim_cell.num_orb):
        prim_cell.add_hopping((1, 0), orb_i, orb_i, 1.0)
        prim_cell.add_hopping((0, 1), orb_i, orb_i, 1.0)
    prim_cell.plot(hop_as_arrows=False, hop_color="gray")


if __name__ == "__main__":
    main()
