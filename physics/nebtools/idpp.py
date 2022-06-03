#! /usr/bin/env python
"""Prepare directories and POSCARs for NEB Calculation."""

import os

import numpy as np
from scipy.optimize import minimize
from ase.io import read, write
from ase.neb import NEB


def norm(x):
    """
    Calculate the norm of vector.

    :param np.ndarray x: incoming vector
    :return: modulus of x
    :rtype: float
    """
    return np.linalg.norm(x).item()


def select_atoms(image, exclude_atoms=None, offset=0):
    """
    Select atoms from an image.

    :param ase.Atoms image: image from which to select atoms
    :param List[int] exclude_atoms: indices of atoms to exclude
    :param int offset: offset of exclude_atoms, use 0 for ase and 1 for vesta
    :return: indices of selected atoms
    :rtype: List[int]
    """
    if exclude_atoms is None:
        exclude_atoms = []
    selected_atoms = [i for i in range(len(image))
                      if i+offset not in exclude_atoms]
    return selected_atoms


def correct_pbc(ref_image, chk_image, **kwargs):
    """
    Correct the coordinates of atoms wrapped back by periodic boundary
    condition to minimize the mismatch.

    :param ase.Atoms ref_image: reference image
    :param ase.Atoms chk_image: image to correct
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


def align_image(ref_image, chk_image, **kwargs):
    """
    Align image to reference to minimize the mismatch.

    :param ase.Atoms ref_image: reference image
    :param ase.Atoms chk_image: image to align
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


def center_image(image, center=0.5, **kwargs):
    """
    Center the image at given position along z-axis.

    :param ase.Atoms image: image to center
    :param float center: position to center the image
    :param kwargs: arguments for 'select_atoms'
    :return: None. Incoming image is modified.
    """
    selected_atoms = select_atoms(image, **kwargs)
    scaled_pos = image.get_scaled_positions()
    old_center = np.sum(scaled_pos[selected_atoms], axis=0)
    old_center /= len(selected_atoms)
    scaled_pos[:, 2] += (center - old_center[2])
    image.set_scaled_positions(scaled_pos)


def interpolate(initial_image, final_image, num_inter_images, **kwargs):
    """
    Generate a path by interpolation between initial and final images.

    :param ase.Atoms initial_image: initial state
    :param ase.Atoms final_image: final state
    :param int num_inter_images: number of intermediate states
    :param kwargs: arguments for 'neb.interpolate'
    :return: path connecting initial and final images
    :rtype: List[ase.Atoms]
    """
    images = [initial_image]
    images += [initial_image.copy() for _ in range(num_inter_images)]
    images.append(final_image)
    neb = NEB(images)
    neb.interpolate(**kwargs)
    return images


def main():
    # File names and number of intermediate states
    poscar_ini = "POSCAR_ini"
    poscar_fin = "POSCAR_fin"
    num_inter_images = 8

    # Normalization parameters
    _correct_pbc = True
    _align_image = True
    _center_image = True
    free_atoms = [9, 31, 32, 33]
    offset = 0
    center_z = 0.2

    # Interpolation parameters
    method = "idpp"
    mic = True

    # Debugging flags
    debug = False

    # Load initial and final images
    image_ini = read(poscar_ini, index=-1)
    image_fin = read(poscar_fin, index=-1)

    # Correct final image
    select_args = {"exclude_atoms": free_atoms, "offset": offset}
    if _correct_pbc:
        correct_pbc(image_ini, image_fin, **select_args)
    if _align_image:
        align_image(image_ini, image_fin, **select_args)

    # Center initial and final images
    if _center_image:
        center_image(image_ini, center_z, **select_args)
        center_image(image_fin, center_z, **select_args)

    # Interpolate
    mep = interpolate(image_ini, image_fin, num_inter_images,
                      method=method, mic=mic)

    # Output
    if not debug:
        vasp_args = {"vasp5": True, "direct": True}
        os.system("rm -rf 0*")
        for i, image in enumerate(mep):
            dir_name = "%02d" % i
            os.mkdir(dir_name)
            write("%s/POSCAR" % dir_name, image, format="vasp", **vasp_args)
    write("mep.traj", mep)


if __name__ == "__main__":
    main()
