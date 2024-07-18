#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import tbplas as tb

from model import make_cell


def run_spin_texture():
    ib = 0
    cell = make_cell(is_qsh=True)

    # Evaluate expectation of sigma_z as scalar field
    k_grid = 2 * (tb.gen_kmesh((240, 240, 1)) - 0.5)
    k_grid[:, 2] = 0.0
    spin_texture = tb.SpinTexture(cell, k_grid, spin_major=False)
    k_cart = spin_texture.k_cart
    sz = spin_texture.eval("z")
    vis = tb.Visualizer()
    vis.plot_scalar(x=k_cart[:, 0], y=k_cart[:, 1], z=sz[:, ib],
                    num_grid=(480, 480), cmap="jet")

    # Evaluate expectations of sigma_x and sigma_y as scalar field
    k_grid = 2 * (tb.gen_kmesh((48, 48, 1)) - 0.5)
    k_grid[:, 2] = 0.0
    spin_texture.k_grid = k_grid
    k_cart = spin_texture.k_cart
    sx = spin_texture.eval("x")
    sy = spin_texture.eval("y")
    vis.plot_vector(x=k_cart[:, 0], y=k_cart[:, 1],
                    u=sx[:, ib], v=sy[:, ib], cmap="jet")

    # Evaluate spin texture of specific energy range
    e_min, e_max = -2, -1.9
    bands = cell.calc_bands(k_grid)[1]
    ik_filt, sx_filt, sy_filt = [], [], []
    for i_k in range(bands.shape[0]):
        idx = np.where((bands[i_k] >= e_min)&(bands[i_k] <= e_max))[0]
        ik_filt.extend([i_k for _ in idx])
        sx_filt.extend(sx[i_k, idx])
        sy_filt.extend(sy[i_k, idx])
    vis.plot_vector(x=k_cart[ik_filt, 0], y=k_cart[ik_filt, 1],
                    u=sx_filt, v=sy_filt, cmap="jet")


def run_fat_band():
    cell = make_cell(is_qsh=True)

    # Evaluate band structure
    k_points = np.array([
        [0.0, 0.0, 0.0],    # Gamma
        [0.5, 0.0, 0.0],    # M
        [2./3, 1./3, 0.0],  # K
        [0.0, 0.0, 0.0],    # Gamma
        [1./3, 2./3, 0.0],  # K'
        [0.5, 0.5, 0.0],    # M'
        [0.0, 0.0, 0.0],    # Gamma
    ])
    k_label = [r"$\Gamma$", "M", "K", r"$\Gamma$",
               r"$K^{\prime}$", r"$M^{\prime}$", r"$\Gamma$"]
    k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40, 40, 40, 40])
    k_len, bands = cell.calc_bands(k_path)

    # Evaluate sigma_? expactation
    component = "y"
    spin_texture = tb.SpinTexture(cell, k_path, spin_major=False)
    projection = spin_texture.eval(component)

    # Plot fat band
    for ib in range(bands.shape[1]):
        plt.scatter(k_len, bands[:, ib], c=projection[:, ib], s=5.0, cmap="jet")
    for x in k_len[k_idx]:
        plt.axvline(x, color="k", linewidth=0.8, linestyle="-")
    for y in (-2.0, -1.9):
        plt.axhline(y, color="k", linewidth=1.0, linestyle=":")
    plt.xlim(k_len[0], k_len[-1])
    plt.xticks(k_len[k_idx], k_label)
    plt.ylabel("Energy (t)")
    plt.colorbar(label=fr"$\langle\sigma_{component}\rangle$")
    plt.show()


if __name__ == "__main__":
    run_spin_texture()
    run_fat_band()
