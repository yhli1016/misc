#! /usr/bin/env python
"""Calculate band structure of TBG using exact diagonalization and TBPM."""

import os
from typing import Tuple

import numpy as np

import tbplas as tb

from tbg import make_tbg


def make_sample(i: int = 5) -> tb.Sample:
    merged_cell = make_tbg(i=i)
    super_cell = tb.SuperCell(merged_cell, dim=(1, 1, 1),
                              pbc=(True, True, False))
    sample = tb.Sample(super_cell)
    print(f"prim_cell num_orb = {merged_cell.num_orb}")
    print(f"prim_cell num_hop = {merged_cell.num_hop}")
    print(f"sample num_orb = {sample.num_orb}")
    print(f"sample num_hop = {sample.num_hop}")
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


def run_diag(enable_mpi: bool = True) -> None:
    sample = make_sample(i=100)
    k_path, k_idx = gen_kpath()
    solver = tb.DiagSolver(sample, enable_mpi=enable_mpi)
    timer = tb.Timer()
    timer.tic("diag")
    k_len, bands = solver.calc_bands(k_path, solver="arpack", which="SM",
                                     k=200)[:2]
    timer.toc("diag")    
    if solver.is_master:
        data_dir = "tbg_diag"
        os.system(f"rm -rf {data_dir} && mkdir {data_dir}")
        np.save(f"{data_dir}/energy", bands)
        np.save(f"{data_dir}/k_len", k_len)
        np.save(f"{data_dir}/k_idx", k_idx)
        timer.report_total_time()


def run_tbpm(enable_mpi: bool = False) -> None:
    sample = make_sample(i=100)

    # Generate k-path along highly-symmetric points
    k_path, k_idx = gen_kpath()
    k_len = tb.gen_kdist(sample.sc0.sc_lat_vec, k_path)

    # Evaluate DOS for each k-point
    config = tb.Config()
    config.generic['nr_time_steps'] = 4096
    config.generic['nr_random_samples'] = 1
    solver = tb.Solver(sample, config, enable_mpi=enable_mpi)
    analyzer = tb.Analyzer(sample, config, enable_mpi=enable_mpi)
    energy_all,  dos_all = [], []
    timer = tb.Timer()
    for i_k, kpt in enumerate(k_path):
        timer.tic("tbpm")
        sample.reset_array()
        sample.set_k_point(kpt)
        sample.rescale_ham()
        corr_dos = solver.calc_corr_dos()
        eng, dos = analyzer.calc_dos(corr_dos)
        energy_all.append(eng)
        dos_all.append(dos)
        timer.toc("tbpm")
    energy_all = np.array(energy_all)
    dos_all = np.array(dos_all)
    if solver.is_master:
        data_dir = "tbg_tbpm"
        os.system(f"rm -rf {data_dir} && mkdir {data_dir}")
        np.save(f"{data_dir}/energy", energy_all)
        np.save(f"{data_dir}/dos", dos_all)
        np.save(f"{data_dir}/k_len", k_len)
        np.save(f"{data_dir}/k_idx", k_idx)
        timer.report_total_time()


if __name__ == "__main__":
    run_diag()
    #run_tbpm()
