"""Functions for NEB calculation."""

import os
from typing import List

import numpy as np
from scipy.optimize import minimize
from ase import Atoms
from ase.io import read, write
from ase.neb import NEB

from .base import get_int


def norm(x: np.ndarray) -> float:
    """
    Calculate the norm of vector.

    :param x: incoming vector
    :return: modulus of x
    """
    return np.linalg.norm(x).item()


def select_atoms(image: Atoms,
                 exclude_atoms: List[int] = None,
                 offset: int = 0) -> List[int]:
    """
    Select atoms from an image.

    :param image: image from which to select atoms
    :param exclude_atoms: indices of atoms to exclude
    :param offset: offset of exclude_atoms, use 0 for ase and 1 for vesta
    :return: indices of selected atoms
    """
    if exclude_atoms is None:
        exclude_atoms = []
    selected_atoms = [i for i in range(len(image))
                      if i+offset not in exclude_atoms]
    return selected_atoms


def correct_pbc(ref_image: Atoms, chk_image: Atoms, **kwargs) -> None:
    """
    Correct the coordinates of atoms wrapped back by periodic boundary
    condition to minimize the mismatch.

    :param ref_image: reference image
    :param chk_image: the image to correct
    :param kwargs: arguments for 'select_atoms'
    :return: None. chk_image is modified.
    """
    # Extract lattice vectors and atomic positions
    ref_cell = ref_image.cell
    ref_pos = ref_image.get_positions()
    chk_pos = chk_image.get_positions()
    diff_pos = ref_pos - chk_pos

    # Select atoms
    selected_atoms = select_atoms(ref_image, **kwargs)

    # Check and fix atoms wrapped back by periodic boundary condition
    neighbor_r0 = [(i, j, k)
                   for i in range(-1, 2)
                   for j in range(-1, 2)
                   for k in range(-1, 2)]
    neighbor_r0 = ref_cell.cartesian_positions(neighbor_r0)
    for i_a in selected_atoms:
        dist = [norm(diff_pos[i_a] - r0) for r0 in neighbor_r0]
        chk_pos[i_a] += neighbor_r0[dist.index(min(dist))]
    chk_image.set_positions(chk_pos)


def align_image(ref_image: Atoms, chk_image: Atoms, **kwargs) -> None:
    """
    Align image to reference to minimize the mismatch.

    :param ref_image: reference image
    :param chk_image: the image to align
    :param kwargs: arguments for 'select_atoms'
    :return: None. chk_image is modified.
    """
    # Extract lattice vectors and atomic positions
    ref_pos = ref_image.get_positions()
    chk_pos = chk_image.get_positions()
    diff_pos = ref_pos - chk_pos

    # Select atoms
    selected_atoms = select_atoms(ref_image, **kwargs)

    # Find the minimal displacement vector
    def _min_obj(x):
        return sum([norm(diff_pos[i_a] - x)**2 for i_a in selected_atoms])

    x0 = np.sum(diff_pos[selected_atoms], axis=0) / len(selected_atoms)
    result = minimize(_min_obj, x0)
    print(result)

    # Align the image according to displacement vector
    if result.status == 0:
        chk_pos += result.x
    else:
        chk_pos += x0
    chk_image.set_positions(chk_pos)


def center_image(image: Atoms, center: float = 0.5, **kwargs) -> None:
    """
    Center the image at given position along z-axis.

    :param image: image to center
    :param center: position to center the image
    :param kwargs: arguments for 'select_atoms'
    :return: None. Incoming image is modified.
    """
    selected_atoms = select_atoms(image, **kwargs)
    scaled_pos = image.get_scaled_positions()
    old_center = np.sum(scaled_pos[selected_atoms], axis=0)
    old_center /= len(selected_atoms)
    scaled_pos[:, 2] += (center - old_center[2])
    image.set_scaled_positions(scaled_pos)


def diff_image(ref_image: Atoms, chk_image: Atoms) -> float:
    """
    Evaluate the difference between given image and reference.

    :param ref_image: reference image
    :param chk_image: the image to diff
    :return: total difference
    """
    ref_pos = ref_image.get_positions()
    chk_pos = chk_image.get_positions()
    diff_pos = ref_pos - chk_pos
    return sum([norm(diff_pos[i_a])**2 for i_a in range(len(ref_image))])


def interpolate(initial_image: Atoms,
                final_image: Atoms,
                num_inter_images: int,
                **kwargs) -> List[Atoms]:
    """
    Generate a path by interpolation between initial and final images.

    :param initial_image: initial state
    :param final_image: final state
    :param num_inter_images: number of intermediate states
    :param kwargs: arguments for 'neb.interpolate'
    :return: path connecting initial and final images
    """
    images = [initial_image]
    images += [initial_image.copy() for _ in range(num_inter_images)]
    images.append(final_image)
    neb = NEB(images)
    neb.interpolate(**kwargs)
    return images


def mep2pos(mep: List[Atoms]) -> None:
    """
    Write mep to POSCAR.

    :param mep: mep to save
    :return: None
    """
    vasp_args = {"vasp5": True, "direct": True}
    os.system("rm -rf 0*")
    for i, image in enumerate(mep):
        dir_name = f"{i:02d}"
        os.mkdir(dir_name)
        write(f"{dir_name}/POSCAR", image, format="vasp", **vasp_args)


def pos2mep() -> None:
    """Convert CONTCAR/POSCAR to mep after NEB calculation."""
    num_image = get_int("INCAR", "IMAGES")
    last_run = get_int("run_neb.sh", "run")
    images = []
    for i in range(num_image+2):
        try:
            image = read(f"{i:02d}/CONTCAR", index=-1)
        except FileNotFoundError:
            image = read(f"{i:02d}/POSCAR", index=-1)
        images.append(image)
    write(f"mep_{last_run}.traj", images)
