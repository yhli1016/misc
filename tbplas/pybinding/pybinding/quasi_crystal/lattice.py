"""Functions for lattice operations."""

from math import sin, cos, sqrt, pi

import numpy as np


INFINITESIMAL = 1.0e-15


__all__ = ["gen_lattice_vectors", "gen_reciprocal_vectors",
           "frac2cart", "cart2frac", "frac2cart_s", "cart2frac_s",
           "rotate_coord", "get_lattice_area", "get_lattice_volume"]


def gen_lattice_vectors(a: float = 1.0,
                        b: float = 1.0,
                        c: float = 1.0,
                        alpha: float = 90.0,
                        beta: float = 90.0,
                        gamma: float = 90.0) -> np.ndarray:
    """
    Generate lattice vectors from given lattice parameters.

    Reference:
    http://www.quantum-espresso.org/Doc/INPUT_PW.html

    :param a: lattice constant 'a' in ANY unit
    :param b: lattice constant 'b' in ANY unit
    :param c: lattice constant 'c' in ANY unit
    :param alpha: angle between a2 and a3 in DEGREE
    :param beta: angle between a3 and a1 in DEGREE
    :param gamma: angle between a1 and a2 in DEGREE
    :return: (3, 3) float64 array
        Cartesian coordinates of lattice vectors in the same unit as a/b/c,
        with each ROW being a vector
    """
    # Check for insane geometry parameters
    for f in (a, b, c, alpha, beta, gamma):
        assert f > INFINITESIMAL

    # Calculate lattice vectors
    alpha = alpha / 180 * pi
    beta = beta / 180 * pi
    gamma = gamma / 180 * pi
    lat_vec = np.zeros((3, 3))
    lat_vec[0, :] = [a, 0, 0]
    lat_vec[1, :] = [b*cos(gamma), b*sin(gamma), 0]
    lat_vec[2, 0] = c * cos(beta)
    lat_vec[2, 1] = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma)
    lat_vec[2, 2] = c * sqrt(1 + 2*cos(alpha)*cos(beta)*cos(gamma)
                             - cos(alpha)**2 - cos(beta)**2
                             - cos(gamma)**2) / sin(gamma)

    # Checks the volume
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    return lat_vec


def gen_reciprocal_vectors(lat_vec: np.ndarray) -> np.ndarray:
    """
    Generate reciprocal lattice vectors from real-space lattice vectors.

    We evaluate reciprocal lattice vectors following
        dot_product(a_i, b_j) = 2 * pi * delta_{ij}
    The formulae based on cross-products are not robust in some cases.

    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of real-space lattice vectors, with each ROW
        being a vector
    :return: (3, 3) float64 array
        Cartesian coordinates of reciprocal lattice vectors, with each ROW
        being a vector. The unit is inverse to the unit of real-space lattice
        vectors.
    """
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    recip_vec = np.zeros((3, 3))
    product = 2 * pi * np.eye(3)
    for i in range(3):
        recip_vec[i] = np.linalg.solve(lat_vec, product[i])
    return recip_vec


def cart2frac(lat_vec: np.ndarray,
              cart_coord: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert Cartesian coordinates to fractional coordinates.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :param cart_coord: (num_coord, 3) float64 array
        Cartesian coordinates to convert
    :param origin: (3, ) float64 array
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        fractional coordinates in basis of lattice vectors
    """
    if not isinstance(lat_vec, np.ndarray):
        lat_vec = np.array(lat_vec)
    if not isinstance(cart_coord, np.ndarray):
        cart_coord = np.array(cart_coord)
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    # Explicit data type is essential, since cart_coord may be an integer matrix
    # (e.g., when making twisted bilayer hetero-structures).
    frac_coord = np.zeros_like(cart_coord, dtype=np.float64)
    conv_mat = np.linalg.inv(lat_vec.T)
    for i, row in enumerate(cart_coord):
        frac_coord[i] = np.matmul(conv_mat, row - origin)
    return frac_coord


def cart2frac_s(lat_vec: np.ndarray,
                cart_coord: np.ndarray,
                origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert a single Cartesian coordinate to fractional.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :param cart_coord: (3, ) float64 array
        Cartesian coordinate to convert
    :param origin: (3, ) float64 array
        Cartesian coordinate of lattice origin
    :return: (3, ) float64 array
        fractional coordinate in basis of lattice vectors
    """
    if not isinstance(lat_vec, np.ndarray):
        lat_vec = np.array(lat_vec)
    if not isinstance(cart_coord, np.ndarray):
        cart_coord = np.array(cart_coord)
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    conv_mat = np.linalg.inv(lat_vec.T)
    frac_coord = np.matmul(conv_mat, cart_coord - origin)
    return frac_coord


def frac2cart(lat_vec: np.ndarray,
              frac_coord: np.ndarray,
              origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert fractional coordinates to Cartesian coordinates.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :param frac_coord: (num_coord, 3) float64 array
        fractional coordinates to convert in basis of lattice vectors
    :param origin: (3, ) float64 array
        Cartesian coordinate of lattice origin
    :return: (num_coord, 3) float64 array
        Cartesian coordinates converted from fractional coordinates
    """
    if not isinstance(lat_vec, np.ndarray):
        lat_vec = np.array(lat_vec)
    if not isinstance(frac_coord, np.ndarray):
        frac_coord = np.ndarray(frac_coord)
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    # Explicit data type is essential, since cart_coord may be an integer matrix
    # (e.g., when making twisted bilayer hetero-structures).
    cart_coord = np.zeros_like(frac_coord, dtype=np.float64)
    conv_mat = np.copy(lat_vec.T, order="C")
    for i, row in enumerate(frac_coord):
        cart_coord[i] = np.matmul(conv_mat, row) + origin
    return cart_coord


def frac2cart_s(lat_vec: np.ndarray,
                frac_coord: np.ndarray,
                origin: np.ndarray = np.zeros(3)) -> np.ndarray:
    """
    Convert a single fractional coordinate to Cartesian.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :param frac_coord: (3, ) float64 array
        fractional coordinate to convert in basis of lattice vectors
    :param origin: (3, ) float64 array
        Cartesian coordinate of lattice origin
    :return: (3, ) float64 array
        Cartesian coordinate converted from fractional coordinate
    """
    if not isinstance(lat_vec, np.ndarray):
        lat_vec = np.array(lat_vec)
    if not isinstance(frac_coord, np.ndarray):
        frac_coord = np.ndarray(frac_coord)
    assert get_lattice_volume(lat_vec) > INFINITESIMAL
    conv_mat = np.copy(lat_vec.T, order="C")
    cart_coord = np.matmul(conv_mat, frac_coord) + origin
    return cart_coord


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
    :param center: (3, ) float64 array
        Cartesian coordinate of rotation center
    :return: (num_coord, 3) float64 array
        rotated Cartesian coordinates
    """
    assert len(center) == 3
    assert axis in ("x", "y", "z")
    if not isinstance(coord, np.ndarray):
        coord = np.array(coord)
    if not isinstance(center, np.ndarray):
        center = np.array(center)

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
    rot_mat = np.copy(rot_mat, order="C")

    # Rotate coordinates
    # Explicit data type is essential, since cart_coord may be an integer matrix
    # (e.g., when making twisted bilayer hetero-structures).
    coord_rot = np.zeros_like(coord, dtype=np.float64)
    for i in range(coord.shape[0]):
        coord_rot[i] = np.matmul(rot_mat, coord[i] - center) + center
    return coord_rot


def get_lattice_area(lat_vec: np.ndarray, direction: str = "c") -> float:
    """
    Calculate the area along given direction.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :param direction: direction along which the area is evaluated, should be in
        ("a", "b", "c")
    :return: area along given direction in squared unit of lattice vectors
    """
    assert direction in ("a", "b", "c")
    if direction == "a":
        a0 = lat_vec[1]
        a1 = lat_vec[2]
    elif direction == "b":
        a0 = lat_vec[2]
        a1 = lat_vec[0]
    else:
        a0 = lat_vec[0]
        a1 = lat_vec[1]
    area = np.linalg.norm(np.cross(a0, a1)).item()
    assert area > INFINITESIMAL
    return area


def get_lattice_volume(lat_vec: np.ndarray) -> float:
    """
    Calculate the volume formed by lattice vectors.
    :param lat_vec: (3, 3) float64 array
        Cartesian coordinates of lattice vectors, with each ROW being a vector
    :return: lattice volume in cubic unit of lattice vectors
    """
    a0 = lat_vec[0]
    a1 = lat_vec[1]
    a2 = lat_vec[2]
    volume = np.abs(np.dot(np.cross(a0, a1), a2)).item()
    assert volume > INFINITESIMAL
    return volume
