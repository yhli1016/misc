#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import tbplas as tb

from model import make_cell


def main(is_qsh: bool = True):
    # Set up primitive cell
    prim_cell = make_cell(is_qsh)

    # Rotate 30 degree counter-clockwise
    tb.spiral_prim_cell(prim_cell, angle=np.pi/6)

    # Reshape to rectangular shape
    lat_sc = np.array([[1, 0, 0], [-1, 2, 0], [0, 0, 1]])
    prim_cell = tb.reshape_prim_cell(prim_cell, lat_sc)

    # Extend along y-direction to increase width
    prim_cell = tb.extend_prim_cell(prim_cell, dim=(1, 3, 1))

    # Disable PBC along y to yield zigzag nano ribbon
    prim_cell.apply_pbc((True, False, False))
    prim_cell.plot()

    # Calculate and plot fat band structure
    k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [0.5, 0.0, 0.0],    # M
        [0.0, 0.0, 0.0],    # M
    ])
    k_label = ["M", r"$\Gamma$", "M"]
    k_path, k_idx = tb.gen_kpath(k_points, [40, 40])
    k_len, bands = prim_cell.calc_bands(k_path)
    component = "z"
    spin_texture = tb.SpinTexture(prim_cell, k_path, spin_major=False)
    projection = spin_texture.eval(component)

    for ib in range(bands.shape[1]):
        plt.scatter(k_len, bands[:, ib], c=projection[:, ib], s=5.0, cmap="jet")
    for x in k_len[k_idx]:
        plt.axvline(x, color="k", linewidth=0.8, linestyle="-")
    plt.xlim(k_len[0], k_len[-1])
    plt.xticks(k_len[k_idx], k_label)
    plt.ylabel("Energy (t)")
    plt.colorbar(label=fr"$\langle\sigma_{component}\rangle$")
    plt.show()


if __name__ == '__main__':
    main()
