import math

import ase
import numpy as np
from ase.io import read, write

import tbplas as tb
import tbmod as tb2


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


def plot_model(atoms: ase.Atoms) -> None:
    """Visualize rotated structure using TBPLaS."""
    prim_cell = tb.PrimitiveCell(atoms.cell)
    print(np.array(atoms.cell))
    for atom in atoms[:2]:
        prim_cell.add_orbital(atom.scaled_position, energy=1.0)
    for atom in atoms[2:]:
        prim_cell.add_orbital(atom.scaled_position, energy=-1.0)
    neighbors = tb.find_neighbors(prim_cell, a_max=1, b_max=1, c_max=1,
                                  max_distance=0.2)
    for term in neighbors:
        prim_cell.add_hopping(term.rn, term.pair[0], term.pair[1], energy=1.0)
    prim_cell.plot(view="ab", with_conj=False, hop_as_arrows=False)
    prim_cell.plot(view="bc", with_conj=False, hop_as_arrows=False)
    prim_cell.plot(view="ca", with_conj=False, hop_as_arrows=False)


def main():
    atoms = read("CONTCAR", index=-1)
    lat_vec = atoms.cell

    # First approach which rotate twice
    angle = -0.5 * math.pi
    axis = np.matmul(lat_vec, [1, 1, 0])
    center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
    lat_vec = rotate_lattice(lat_vec, angle, axis, center)
    angle = -0.25 * math.pi
    axis = np.array([0, 0, 1], dtype=np.float64)
    center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
    lat_vec = rotate_lattice(lat_vec, angle, axis, center)

    # # Second approach which rotates only once
    # angle = -0.25 * math.pi
    # axis = np.matmul(lat_vec, [-1, 1, 0])
    # center = np.matmul(lat_vec, [0.5, 0.5, 0.5])
    # lat_vec = rotate_lattice(lat_vec, angle, axis, center)

    # Save structure
    atoms.set_cell(lat_vec, scale_atoms=True)
    write("POSCAR", atoms)
    plot_model(atoms)


if __name__ == "__main__":
    main()
