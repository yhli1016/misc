import math

import numpy as np

import tbplas as tb


def make_prim_cell(origin: np.ndarray = np.zeros(3)) -> tb.PrimitiveCell:
    prim_cell = tb.PrimitiveCell(lat_vec=np.eye(3), origin=origin)
    prim_cell.add_orbital((0, 0), 0.0)
    prim_cell.add_hopping((1, 0), 0, 0, 1.0)
    prim_cell.add_hopping((0, 1), 0, 0, 1.0)
    return prim_cell


def main():
    scatter_cell = make_prim_cell()
    scatter = tb.SuperCell(scatter_cell, (10, 10, 1), (False, False, False))

    left_lead_cell = make_prim_cell(origin=(-3.0, 0.0, 0.0))
    left_lead = tb.SuperCell(left_lead_cell, (3, 10, 1), (False, False, False))

    right_lead_cell = make_prim_cell(origin=(10.0, 0.0, 0.0))
    right_lead = tb.SuperCell(right_lead_cell, (3, 10, 1), (False, False, False))

    # Build inter-cell hopping terms
    inter_hop0 = tb.SCInterHopping(scatter, left_lead)
    neighbors = tb.find_neighbors(scatter, left_lead, a_max=0, b_max=0, max_distance=0.11)
    for term in neighbors:
        i, j = term.pair
        inter_hop0.add_hopping(term.rn, i, j, -0.5)

    inter_hop1 = tb.SCInterHopping(scatter, right_lead)
    neighbors = tb.find_neighbors(scatter, right_lead, a_max=0, b_max=0, max_distance=0.11)
    for term in neighbors:
        i, j = term.pair
        inter_hop1.add_hopping(term.rn, i, j, -0.5)

    sample = tb.Sample(scatter, left_lead, right_lead, inter_hop0, inter_hop1)
    sample.plot(with_orbitals=True, with_conj=False, hop_as_arrows=False, with_cells=False,
                sc_hop_colors=["g", "r", "b"], inter_hop_colors=["k", "k"])


if __name__ == "__main__":
    main()
