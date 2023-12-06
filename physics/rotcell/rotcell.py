#! /usr/bin/env python
import math

import ase
import numpy as np
from numpy.linalg import norm
from ase.io import read, write


def rotate_coord(coord: np.ndarray,
                 angle: float = 0.0,
                 axis: np.ndarray = np.array([0, 0, 1]),
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
    :param axis: vector
        Cartesian coordinate of the rotation axis
    :param center: vector
        Cartesian coordinate of rotation center
    :return: (num_coord, 3) float64 array
        rotated Cartesian coordinates
    """
    if not isinstance(coord, np.ndarray):
        coord = np.array(coord)
    if not isinstance(center, np.ndarray):
        center = np.array(center)
    if len(center) != 3:
        raise ValueError(f"Length of rotation center should be 3")

    # Determine rotation matrix
    axis /= np.linalg.norm(axis)
    ux, uy, uz = axis
    u_prod = np.array([[0, -uz, uy], [uz, 0, -ux], [-uy, ux, 0]])
    u_tens = np.tensordot(axis, axis, axes=0)
    cos_ang, sin_ang = math.cos(angle), math.sin(angle)
    rot_mat = cos_ang * np.eye(3) + sin_ang * u_prod + (1 - cos_ang) * u_tens

    # Rotate coordinates
    coord_rot = np.zeros(shape=coord.shape, dtype=coord.dtype)
    for i in range(coord.shape[0]):
        coord_rot[i] = np.matmul(rot_mat, coord[i] - center) + center
    return coord_rot


def rotate_lattice(lattice: np.ndarray,
                   angle: float,
                   axis: np.ndarray,
                   center: np.ndarray) -> np.ndarray:
    """Rotate lattice vectors according to Euler angles."""
    end_points = np.vstack((np.zeros(3), lattice))
    end_points = rotate_coord(end_points, angle=angle, center=center, axis=axis)
    lattice = end_points[1:] - end_points[0]
    return lattice


def calc_norm_vector(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Calculate the norm vector of the plane formed by v1 and v2."""
    v3 = np.cross(v1, v2)
    v3 /= norm(v3)
    return v3


def calc_angle(v1: np.ndarray, v2: np.ndarray) -> float:
    """Calculate the angle between v2 and v2."""
    angle = math.acos(np.dot(v1, v2) / norm(v1) / norm(v2))
    return angle


def rotate_atoms(atoms: ase.Atoms,
                 idx0: int,
                 idx1: int,
                 ref_vec: np.ndarray) -> None:
    """
    Rotate the atoms such that the given bond parallels to reference vector.

    :param atoms: atoms to rotate
    :param idx0: starting atom index of the bond
    :param idx1: ending atom index of the bond
    :param ref_vec: Cartesian coordinates of the reference vector
    :return: None. The incoming image is modified.
    """
    lat_vec = atoms.cell
    atom_pos = atoms.get_positions()
    ra, rb = atom_pos[idx0], atom_pos[idx1]
    ab = rb - ra
    angle = calc_angle(ab, ref_vec)
    axis = calc_norm_vector(ab, ref_vec)
    center = ra
    lat_vec = rotate_lattice(lat_vec, angle, axis, center)
    atoms.set_cell(lat_vec, scale_atoms=True)


def test_r():
    atoms = read("POSCAR_R.origin", index=-1)
    algo = 0

    if algo == 0:
        # First approach which rotate twice
        lat_vec = atoms.cell
        angle = -0.5 * math.pi
        axis = np.matmul(lat_vec, [1, 1, 0])
        center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
        lat_vec = rotate_lattice(lat_vec, angle, axis, center)
        angle = -0.25 * math.pi
        axis = np.array([0, 0, 1], dtype=np.float64)
        center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
        lat_vec = rotate_lattice(lat_vec, angle, axis, center)
        atoms.set_cell(lat_vec, scale_atoms=True)

    elif algo == 1:
        # Alternative of first approach
        idx0, idx1 = 1, 5
        ref_vec = np.array([0.0, 0.0, 1.0])
        rotate_atoms(atoms, idx0, idx1, ref_vec)
        idx0, idx1 = 1, 2
        ref_vec = np.array([-1.0, -1.0, 0.0])
        rotate_atoms(atoms, idx0, idx1, ref_vec)

    elif algo == 2:
        # Second approach which rotates only once
        lat_vec = atoms.cell
        angle = -0.25 * math.pi
        axis = np.matmul(lat_vec, [-1, 1, 0])
        center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
        lat_vec = rotate_lattice(lat_vec, angle, axis, center)
        atoms.set_cell(lat_vec, scale_atoms=True)

    else:
        # Alternative of first approach
        idx0, idx1 = 1, 2
        ref_vec = np.array([0.0, 0.0, -1.0])
        rotate_atoms(atoms, idx0, idx1, ref_vec)

    # Save structure
    write("POSCAR_R.rot", atoms)


def test_m1():
    atoms = read("POSCAR_M1.origin", index=-1)
    idx0, idx1 = 2, 9
    ref_vec = np.array([0.0, 0.0, 1.0])
    rotate_atoms(atoms, idx0, idx1, ref_vec)
    write("POSCAR_M1.rot", atoms)


if __name__ == "__main__":
    test_r()
    test_m1()
