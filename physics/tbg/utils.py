from math import cos, sin, floor, ceil

import numpy as np
from ase import Atoms


def cart2frac(lattice_vectors: np.ndarray,
              cartesian_coordinates: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert Cartesian coordinates to fractional coordinates.

    :param lattice_vectors: (3, 3) float64 array
        Cartesian coordinates of lattice vectors
    :param cartesian_coordinates: (num_coord, 3) float64 array
        Cartesian coordinates to convert
    :param origin: float64 array of length 3
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        fractional coordinates in basis of lattice vectors
    """
    if not isinstance(lattice_vectors, np.ndarray):
        lattice_vectors = np.array(lattice_vectors)
    if not isinstance(cartesian_coordinates, np.ndarray):
        cartesian_coordinates = np.array(cartesian_coordinates)
    fractional_coordinates = np.zeros(cartesian_coordinates.shape)
    conversion_matrix = np.linalg.inv(lattice_vectors.T)
    for i, row in enumerate(cartesian_coordinates):
        fractional_coordinates[i] = np.matmul(conversion_matrix,
                                              (row - origin).T)
    return fractional_coordinates


def frac2cart(lattice_vectors: np.ndarray,
              fractional_coordinates: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert fractional coordinates to Cartesian coordinates.

    :param lattice_vectors: (3, 3) float64 array
        Cartesian coordinates of lattice vectors
    :param fractional_coordinates: (num_coord, 3) float64 array
        fractional coordinates to convert in basis of lattice vectors
    :param origin: float64 array of length 3
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        Cartesian coordinates converted from fractional coordinates
    """
    if not isinstance(lattice_vectors, np.ndarray):
        lattice_vectors = np.array(lattice_vectors)
    if not isinstance(fractional_coordinates, np.ndarray):
        fractional_coordinates = np.ndarray(fractional_coordinates)
    cartesian_coordinates = np.zeros(fractional_coordinates.shape)
    conversion_matrix = lattice_vectors.T

    for i, row in enumerate(fractional_coordinates):
        cartesian_coordinates[i] = np.matmul(conversion_matrix, row.T) + origin
    return cartesian_coordinates


def rotate_coord(coord: np.ndarray,
                 angle: float = 0.0,
                 axis: str = "z",
                 center: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Rotate Cartesian coordinates according to Euler angles.

    Reference:
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    Note that the right-hand rule is used, i.e., a positive value means the
    angle is counter-clockwise.

    :param coord: (num_coord, 3) float64 array
        Cartesian coordinates to rotate
    :param angle: float
        rotation angle in RADIAN, not degrees
    :param axis: string
        axis around which the rotation is performed
    :param center: float64 vector
        Cartesian coordinate of rotation center
    :return: (num_coord, 3) float64 array
        rotated Cartesian coordinates
    :raises ValueError: if axis is not "x", "y" or "z"
    """
    if not isinstance(coord, np.ndarray):
        coord = np.array(coord)
    if not isinstance(center, np.ndarray):
        center = np.array(center)
    if len(center) != 3:
        raise ValueError(f"Length of rotation center should be 3")
    if axis not in ("x", "y", "z"):
        raise ValueError("Axis should be in 'x', 'y', 'z'")

    # Determine rotation matrix
    if axis == "x":
        u = np.array([1, 0, 0])
    elif axis == "y":
        u = np.array([0, 1, 0])
    else:
        u = np.array([0, 0, 1])
    ux, uy, uz = u
    u_prod = np.array([[0, -uz, uy], [uz, 0, -ux], [-uy, ux, 0]])
    u_tens = np.tensordot(u, u, axes=0)
    cos_ang, sin_ang = cos(angle), sin(angle)
    rot_mat = cos_ang * np.eye(3) + sin_ang * u_prod + (1 - cos_ang) * u_tens

    # Rotate coordinates
    coord_rot = np.zeros(shape=coord.shape, dtype=coord.dtype)
    for i in range(coord.shape[0]):
        coord_rot[i] = np.matmul(rot_mat, coord[i] - center) + center
    return coord_rot


def twist_cell(prim_cell: Atoms,
               angle: float = 0.0,
               center: np.ndarray = np.zeros(3),
               shift: float = 0.0) -> None:
    """
    Rotate and shift primitive cell with respect to z-axis.

    NOTE: this function returns nothing. But the incoming primitive cell
    will be modified in-place.

    :param prim_cell: primitive cell to twist
    :param angle: twisting angle in RADIANs, NOT degrees
    :param center: float64 vector
        Cartesian coordinates of the rotation center in Angstrom
    :param shift: distance of shift in Angstrom
    :return: None
    """
    # Get rotated lattice vectors and origin
    end_points = np.vstack((np.zeros(3), prim_cell.cell))
    end_points = rotate_coord(end_points, angle=angle, center=center)
    lat_vec = end_points[1:] - end_points[0]

    # Get rotated atomic positions
    positions = prim_cell.get_positions()
    positions = rotate_coord(positions, angle=angle, center=center)

    # Reset lattice vectors and atomic positions
    prim_cell.set_cell(lat_vec, scale_atoms=False)
    prim_cell.positions = positions + np.array([0, 0, shift])


def reshape_cell(prim_cell: Atoms,
                 lat_frac: np.ndarray,
                 delta: float = 1e-2) -> Atoms:
    """
    Reshape primitive cell to given lattice vectors.

    :param prim_cell: primitive cell from which the reshaped cell is constructed
    :param lat_frac: (3, 3) float64 array
        FRACTIONAL coordinates of lattice vectors of reshaped cell in basis
        vectors of primitive cell
    :param delta: small displacement added to orbital positions such that
        orbitals fall on cell borders will not be clipped
    :return: reshaped supercell
    """
    # Reshaped cell lattice vectors
    lat_cart = np.zeros((3, 3), dtype=np.float64)
    for i_dim in range(3):
        lat_cart[i_dim] = np.matmul(lat_frac[i_dim], prim_cell.cell)

    # Determine searching range
    # sum_vec is actually a0+a1, a1+a2 or a2+a0 depending on j, i.e., the
    # diagonal vector. If it is not taken into consideration, orbitals and
    # hopping terms in the top right corner of reshaped lattice may be missing.
    rn_range = np.zeros((3, 2), dtype=np.int32)
    for i in range(3):
        sum_vec = lat_frac.sum(axis=0) - lat_frac[i]
        for j in range(3):
            rn_range[j, 0] = min(rn_range[j, 0], floor(sum_vec[j]))
            rn_range[j, 1] = max(rn_range[j, 1], ceil(sum_vec[j]))
    rn_range[:, 0] -= 1
    rn_range[:, 1] += 1

    # Conversion matrix of fractional coordinates from primitive cell to
    # reshaped cell: x_res = x_prim * conv_mat, with x_new and x_prim
    # being ROW vectors
    conv_mat = np.linalg.inv(lat_frac)

    # Function for getting cell index from fractional coordinates
    def _get_cell_index(x):
        return floor(x.item(0)), floor(x.item(1)), floor(x.item(2))

    # Add orbitals
    symbols, scaled_positions = [], []
    orb_pos = prim_cell.get_scaled_positions()
    ra_min, ra_max = rn_range.item(0, 0), rn_range.item(0, 1)
    rb_min, rb_max = rn_range.item(1, 0), rn_range.item(1, 1)
    rc_min, rc_max = rn_range.item(2, 0), rn_range.item(2, 1)
    for i_a in range(ra_min, ra_max+1):
        for i_b in range(rb_min, rb_max+1):
            for i_c in range(rc_min, rc_max+1):
                rn = (i_a, i_b, i_c)
                for i_o, pos in enumerate(orb_pos):
                    res_pos = np.matmul(rn + pos, conv_mat)
                    res_rn = _get_cell_index(res_pos + delta)
                    if res_rn == (0, 0, 0):
                        symbols.append(prim_cell.symbols[i_o])
                        scaled_positions.append(res_pos)

    # Create reshaped cell
    res_cell = Atoms(cell=lat_cart, symbols=symbols,
                      scaled_positions=scaled_positions)
    return res_cell


def make_layer(prim_cell: Atoms, hetero_lattice: np.ndarray) -> Atoms:
    """
    Make one layer in the hetero-structure by reshaping primitive cell to
    given lattice vectors.

    :param prim_cell: incoming primitive cell
    :param hetero_lattice: (3, 3) float64 array
        Cartesian coordinates of hetero-structure lattice vectors
    :return: layer in the hetero-structure
    """
    hetero_lattice_frac = cart2frac(prim_cell.cell, hetero_lattice)
    hetero_layer = reshape_cell(prim_cell, hetero_lattice_frac)
    return hetero_layer
