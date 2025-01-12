#! /usr/bin/env python

import numpy as np

import tbplas as tb

from model import make_cell


def main(is_qsh: bool = True):
    # Set up primitive cell
    prim_cell = make_cell(is_qsh)

    # Get evolution of Wannier centers
    ka_array = np.linspace(-0.5, 0.5, 200)
    kb_array = np.linspace(0.0, 0.5, 200)
    kc = 0.0
    z2 = tb.Z2(prim_cell, num_occ=2)
    kb_array, phases = z2.calc_phases(ka_array, kb_array, kc)

    # Visualize the evolution
    vis = tb.Visualizer()
    vis.plot_phases(kb_array, phases, scatter=True)
    vis.plot_phases(kb_array, phases, scatter=True, polar=True)


if __name__ == '__main__':
    main()
