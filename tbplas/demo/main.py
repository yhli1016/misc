#! /usr/bin/env python
import json
import math
import sys
import numpy as np
import tbplas as tb

import fractal
import qs
import tbg


def make_fractal(data: dict) -> tb.Sample:
    lattice = np.eye(3, dtype=np.float64)
    prim_cell = tb.PrimitiveCell(lattice)
    prim_cell.add_orbital((0, 0))
    prim_cell.add_hopping((1, 0), 0, 0, 1.0)
    prim_cell.add_hopping((0, 1), 0, 0, 1.0)
    try:
        fractal_kind = data["fractal_kind"]
    except KeyError:
        fractal_kind = [2, 3, 3]
    model = fractal.top_down(prim_cell, start_width=fractal_kind[0], extension=fractal_kind[1], iteration=fractal_kind[2])
    model = tb.Sample(tb.SuperCell(model, dim=(1, 1, 1), pbc=(False, False, False)))
    return model


def make_quasi_crystal(data: dict) -> tb.Sample:
    try:
        radius = data["quasi_crystal_kind"]
    except KeyError:
        radius = 3.0
    dim_a = int(radius / 0.142 * 1.5)
    dim = (dim_a, dim_a, 1)
    angle = 30 / 180 * math.pi
    center = np.array((2./3, 2./3, 0))
    prim_cell = tb.make_graphene_diamond()
    model = qs.make_quasi_crystal_pc(prim_cell, dim, angle, center, radius)
    model = tb.Sample(tb.SuperCell(model, dim=(1, 1, 1), pbc=(False, False, False)))
    return model


def make_tbg(data: dict) -> tb.Sample:
    try:
        tbg_kind = data["tbg_kind"]
    except KeyError:
        tbg_kind = 5
    model = tb.Sample(tb.SuperCell(tbg.make_tbg(tbg_kind), dim=(3, 3, 1), pbc=(True, True, False)))
    return model


def build(data: dict) -> None:
    # Build model
    if data["model"] == "fractal":
        model = make_fractal(data)
    elif data["model"] == "quasi_crystal":
        model = make_quasi_crystal(data)
    else:
        model = make_tbg(data)

    # Plot
    model.plot(with_cells=False, with_orbitals=False, hop_as_arrows=False,
               with_conj=False, hop_eng_cutoff=0.5, fig_name=data["output"], fig_dpi=100)


def calc_dos(sample: tb.Sample, data) -> None:
    config = tb.Config()
    config.generic['nr_time_steps'] = 1024
    config.generic['nr_random_samples'] = 1
    solver = tb.Solver(sample, config)
    analyzer = tb.Analyzer(sample, config)
    corr_dos = solver.calc_corr_dos()
    eng, dos = analyzer.calc_dos(corr_dos)
    vis = tb.Visualizer()
    vis.plot_dos(eng, dos, fig_name=data["output"])


def calc_wfc(sample: tb.Sample, data:dict) -> None:
    config = tb.Config()
    config.generic['nr_random_samples'] = 1
    config.generic['nr_time_steps'] = 1024
    config.quasi_eigenstates['energies'] = [data["wfc_energy"]]

    solver = tb.Solver(sample, config)
    qs = solver.calc_quasi_eigenstates()
    wfc = np.abs(qs[0])**2
    wfc /= wfc.max()
    vis = tb.Visualizer()
    vis.plot_wfc(sample, wfc, scatter=True, site_size=wfc*100, site_color="cmap", cmap="jet",
                 with_model=False, model_style={"alpha": 0.002, "color": "gray"},
                 fig_name=data["output"], fig_size=(6.4, 4.8), fig_dpi=100)


def calculate(data: dict) -> None:
    # Build model
    if data["model"] == "fractal":
        model = make_fractal(data)
    elif data["model"] == "quasi_crystal":
        model = make_quasi_crystal(data)
    else:
        model = make_tbg(data)
    model.rescale_ham()

    # Do the calculation
    if data["job"] == "dos":
        calc_dos(model, data)
    else:
        calc_wfc(model, data)


def main(args):
    if len(args) > 1:
        path = args[1]
    else:
        path = ""
    with open(path + 'input.json') as f:
        data = json.load(f)

    # Check input
    if data["job"] not in ("build", "dos", "wfc"):
        raise ValueError(f"Unkown job type {data['job']}")
    if data["model"] not in ("fractal", "quasi_crystal", "tbg"):
        raise ValueError(f"Unkown model type {data['model']}")

    # Process the request
    if data["job"] == "build":
        build(data)
    else:
        calculate(data)


if __name__ == "__main__":
    main(sys.argv)
