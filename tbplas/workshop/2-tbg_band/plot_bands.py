#! /usr/bin/env python

import sys

import numpy as np

from plotutils import plot_diag, plot_tbpm


def main(data_dir: str = "demo_diag") -> None:
    if data_dir == "demo_diag":
        e_min, e_max = -8, 8
        num_grid = (1, 1)
        tbpm = False

    elif data_dir == "demo_tbpm":
        e_min, e_max = -10, 10
        num_grid = (1000, 1000)
        tbpm = True

    elif data_dir in ("tbg31_diag/prim_cell_dense", "tbg31_diag/prim_cell_sparse",
                      "tbg31_diag/sample_dense", "tbg31_diag/sample_sparse",
                      "tbg100_diag/sample_sparse"):
        e_min, e_max = -0.21, 0.21
        num_grid = (1, 1)
        tbpm = False

    elif data_dir in ("tbg31_tbpm", "tbg100_tbpm"):
        e_min, e_max = -0.21, 0.21
        num_grid = (4000, 4000)
        tbpm = True

    else:
        raise ValueError(f"Unknown data_dir {data_dir}")

    # Load data
    energy = np.load(f"{data_dir}/energy.npy")
    k_len = np.load(f"{data_dir}/k_len.npy")
    k_idx = np.load(f"{data_dir}/k_idx.npy")
    k_label = ["K", r"$\Gamma$", "M", "K"]

    # Plot
    if tbpm:
        dos = np.load(f"{data_dir}/dos.npy")
        plot_tbpm(energy, dos, k_idx=k_idx, k_len=k_len, k_label=k_label,
                  e_min=e_min, e_max=e_max, num_grid=num_grid, cmap="jet")
    else:
        plot_diag(energy, k_idx=k_idx, k_len=k_len, k_label=k_label,
                  e_min=e_min, e_max=e_max, scatter=True)


if __name__ == "__main__":
    main(sys.argv[1])
