#! /usr/bin/env python
import sys

import numpy as np

import pybinding as pb

from mask import Box, Mask
from utils import Timer


FAIR_PLAY = True


def top_down_polygon(
        lattice: pb.Lattice,
        start_width: int,
        iteration: int, 
        extension: int) -> pb.Model:
    """
    Making fractal using Polygon tool. Fails for iteration > 3.
    """
    # Create the mask
    final_width = start_width * extension**iteration
    start_box = Box(0, 0, final_width - 1, final_width - 1)
    mask = Mask(start_box, num_grid=extension, num_iter=iteration)

    # Create the model
    offset = final_width * 0.5
    full_rect = pb.rectangle(x=final_width, y=final_width)
    for box in mask.boxes:
        if box.void:
            x0, y0 = box.i0, box.j0
            x1, y1 = box.i1, box.j1
            vac = np.array([[x0, y0], [x0, y1], [x1, y1], [x1, y0]]) - offset
            full_rect -= pb.Polygon(vac)
    return pb.Model(lattice, full_rect, pb.translational_symmetry(a1=False, a2=False))


def vacancy(
        final_width: int,
        iteration: int,
        extension: int) -> pb.site_position_modifier:
    """Make site_state_modifier to remove orbitals."""
    offset = final_width * 0.5
    start_box = Box(0, 0, final_width - 1, final_width - 1)
    mask = Mask(start_box, num_grid=extension, num_iter=iteration)

    @pb.site_state_modifier
    def modifier(state, x, y):
        for i in range(x.shape[0]):
            for box in mask.boxes:
                if (box.void and
                    box.i0 <= x[i] + offset - 0.1 < box.i1 and
                    box.j0 <= y[i] + offset < box.j1):
                    state[i] = False
        return state

    return modifier


def top_down_modifier(
        lattice: pb.Lattice,
        start_width: int,
        iteration: int,
        extension: int) -> pb.Model:
    """
    Making fractal using site_state_modifier.
    """
    # Create the mask
    final_width = start_width * extension**iteration
    full_rect = pb.rectangle(x=final_width, y=final_width)
    return pb.Model(lattice, 
                    full_rect,
                    vacancy(final_width, iteration, extension),
                    pb.translational_symmetry(a1=False, a2=False))


def main():
    d = 1
    t = 1
    on_site = 0.1 if FAIR_PLAY else 0.0
    lattice = pb.Lattice(a1=[d, 0], a2=[0, d])
    lattice.add_sublattices(('A', [0, 0], on_site))
    lattice.add_hoppings(
        ([0, 1], 'A', 'A', t),
        ([1, 0], 'A', 'A', t),
        ([1, 1], 'A', 'A', t),
        ([1, -1], 'A', 'A', t),
    )

    model = top_down_modifier(lattice, 2, 3, 3)
    if FAIR_PLAY:
        ham = model.hamiltonian
        print(ham.shape[0], ham.nnz)


if __name__ == "__main__":
    timer = Timer()
    for i in range(5):
        timer.tic("fractal")
        main()
        timer.toc("fractal")
        timer.report_time()
