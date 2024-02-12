#! /usr/bin/env python
"""Prepare directories and POSCARs for NEB Calculation."""

from ase.io import read, write

from nebtools.neb import correct_pbc, align_image, center_image, interpolate, mep2pos


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
    apply_constraint = True

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

    # Interpolate
    mep = interpolate(image_ini, image_fin, num_inter_images,
                      method=method, mic=mic, apply_constraint=apply_constraint)

    # Center images
    if _center_image:
        for image in mep:
            center_image(image, center_z, **select_args)

    # Output
    if not debug:
        mep2pos(mep)
    write("mep_0.traj", mep)


if __name__ == "__main__":
    main()
