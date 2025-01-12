import math
import numpy as np
import tbplas as tb


def make_cell(is_qsh: bool = True) -> tb.PrimitiveCell:
    """
    Set up the primitive cell of Kane-Mele model.

    Ref: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.146802

    :param is_qsh: whether is the model is in quantum spin Hall phase
    :return: the primitive cell of Kane-Mele model
    """
    # Parameters
    lat = 1.
    t = -1.
    lamb_so = 0.06 * t
    lamb_r = 0.05 * t
    if is_qsh:
        lamb_nu = 0.1 * t  # QSH phase
    else:
        lamb_nu = 0.4 * t  # normal insulating phase

    # Function for generating unit vector in xOy plane
    def _ev(_theta):
        _theta *= (math.pi / 180)
        return np.array([math.cos(_theta), math.sin(_theta), 0.0])

    # Lattice vectors
    vectors = np.array([_ev(-30)*lat, _ev(30)*lat, [0, 0, 1]])
    prim_cell = tb.PrimitiveCell(vectors, unit=1.0)

    # Add orbitals
    prim_cell.add_orbital((1./3, 1./3), energy=lamb_nu, label="A+")
    prim_cell.add_orbital((1./3, 1./3), energy=lamb_nu, label="A-")
    prim_cell.add_orbital((2./3, 2./3), energy=-lamb_nu, label="B+")
    prim_cell.add_orbital((2./3, 2./3), energy=-lamb_nu, label="B-")

    # Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_y = np.array([[0, -1j], [1j, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])

    # Add 1st nearest hopping terms, which are ordinary hopping + Rashba SOC
    # dij_table is from sub-lattice A to B.
    dij_table = {(-1, 0, 0): _ev(120), (0, -1, 0): _ev(-120), (0, 0, 0): _ev(0)}
    for rn, dij in dij_table.items():
        hop_mat = t * np.eye(2, dtype=np.complex128)
        hop_mat += 1j * lamb_r * (dij[1] * sigma_x - dij[0] * sigma_y)
        for ii in range(2):
            for jj in range(2):
                prim_cell.add_hopping(rn, ii, jj+2, hop_mat.item(ii, jj))

    # Add 2nd nearest hopping terms, which are inter-site SOC
    # vij_table is for sub-lattice A. For B it is the opposite.
    vij_table = {(-1, 1, 0): -1, (-1, 0, 0): 1, (0, -1, 0): -1}
    for rn, vij in vij_table.items():
        hop_mat = 1j * lamb_so * vij * sigma_z
        for ii in range(2):
            for jj in range(2):
                hij = hop_mat.item(ii, jj)
                prim_cell.add_hopping(rn, ii, jj, hij)
                prim_cell.add_hopping(rn, ii+2, jj+2, -hij)
    return prim_cell
