#! /usr/bin/env python
"""Demonstrate the band structure from exact-diagonalization and TBPM."""

import os
from typing import Tuple

import numpy as np

import tbplas as tb


def make_graphene() -> tb.Sample:
    prim_cell = tb.make_graphene_diamond()
    super_cell = tb.SuperCell(prim_cell, dim=(3, 3, 1),
                              pbc=(True, True, False))
    sample = tb.Sample(super_cell)
    return sample


def gen_kpath() -> Tuple[np.ndarray, np.ndarray]:
    k_points = np.array([
        [2. / 3, 1. / 3, 0.0],
        [0.0, 0.0, 0.0],
        [1. / 2, 0.0, 0.0],
        [2. / 3, 1. / 3, 0.0],
    ])
    k_path, k_idx = tb.gen_kpath(k_points, [40, 40, 40])
    return k_path, k_idx


def run_diag() -> None:
    sample = make_graphene()
    k_path, k_idx = gen_kpath()
    k_len, bands = sample.calc_bands(k_path)
    data_dir = "demo_diag"
    os.system(f"rm -rf {data_dir} && mkdir {data_dir}")
    np.save(f"{data_dir}/energy", bands)
    np.save(f"{data_dir}/k_len", k_len)
    np.save(f"{data_dir}/k_idx", k_idx)


def run_tbpm() -> None:
    sample = make_graphene()

    # Generate k-path along highly-symmetric points
    k_path, k_idx = gen_kpath()
    k_len = tb.gen_kdist(sample.sc0.sc_lat_vec, k_path)

    # Evaluate DOS for each k-point
    config = tb.Config()
    config.generic['nr_time_steps'] = 1024
    config.generic['nr_random_samples'] = 1
    solver = tb.Solver(sample, config)
    analyzer = tb.Analyzer(sample, config)
    energy_all,  dos_all = [], []
    for i_k, kpt in enumerate(k_path):
        sample.reset_array()
        sample.set_k_point(kpt)
        sample.rescale_ham()
        corr_dos = solver.calc_corr_dos()
        eng, dos = analyzer.calc_dos(corr_dos)
        energy_all.append(eng)
        dos_all.append(dos)
    energy_all = np.array(energy_all)
    dos_all = np.array(dos_all)
    if solver.is_master:
        data_dir = "demo_tbpm"
        os.system(f"rm -rf {data_dir} && mkdir {data_dir}")
        np.save(f"{data_dir}/energy", energy_all)
        np.save(f"{data_dir}/dos", dos_all)
        np.save(f"{data_dir}/k_len", k_len)
        np.save(f"{data_dir}/k_idx", k_idx)


if __name__ == "__main__":
    run_diag()
    run_tbpm()
