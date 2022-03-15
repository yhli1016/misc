#! /usr/bin/env python
"""Plot mep to figures via POV-Ray."""

import os
from configparser import ConfigParser

from ase.io import read, write

import idpp


def adjust_ini(prefix, figure_width=800):
    """
    Adjust settings in ini file.

    :param str prefix: prefix of ini file
    :param int figure_width: width of figure
    :return: None
    """
    # Read ini file
    config = ConfigParser()
    with open(f"{prefix}.ini", "r") as ini_file:
        ini_content = ini_file.readlines()
    ini_content.insert(0, "[default]\n")
    config.read_string("".join(ini_content))

    # Adjust settings
    default_width = config.getfloat("default", "Width")
    default_height = config.getfloat("default", "Height")
    aspect_ratio = default_height / default_width
    image_height = int(figure_width * aspect_ratio)
    config.set("default", "Width", str(figure_width))
    config.set("default", "Height", str(image_height))

    # Update ini file
    with open(f"{prefix}.ini", "w") as ini_file:
        config.write(ini_file)


def render(prefix, image, image_width, proj_set, pov_set):
    """
    Render atoms.

    :param str prefix: prefix prepended to output files
    :param ase.Atoms image: image to render
    :param int image_width: width of image in pixel
    :param dict proj_set: general projection settings
    :param dict pov_set: POV-Ray settings
    :return: None
    """
    # Write ini and pov file
    write(f"{prefix}.pov", image, **proj_set, povray_settings=pov_set)
    adjust_ini(prefix, image_width)

    # Run POV-Ray
    # os.system(f"pvengine64 /exit {prefix}[default]")
    os.system(f"povray {prefix}[default]")

    # Clean up
    for suffix in ("ini", "pov"):
        os.remove(f"{prefix}.{suffix}")


def main():
    # TODO: implement 'center_image'
    correct_pbc = True
    align_image = True
    rotate = False
    selected_images = "all"
    style = "pale"  # ase2 ase3 glass simple pale intermediate vmd jmol
    figure_width = 800

    # Load and select images
    images = read("mep.traj", index=":")
    ref_image = images[0]
    if selected_images == "all":
        selected_images = range(len(images))
    images = [images[_] for _ in selected_images]

    # Documentation of the parameters:
    # https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=pov
    proj_set = {
        'show_unit_cell': 0,
    }
    pov_set = {
        'transparent': True,
        'background': 'White',
        'textures': [style for _ in range(len(ref_image))]
    }

    # Correct pbc and align images
    for image in images[1:]:
        if correct_pbc:
            idpp.correct_pbc(ref_image, image, selected_atoms="all")
        if align_image:
            idpp.align_image(ref_image, image, selected_atoms="all")

    # Render images
    for i, image in enumerate(images):
        if rotate:
            # proj_set['rotation'] = "-90x"
            image.rotate(-90, "x", rotate_cell=True)
        render(f"{selected_images[i]}", image, figure_width, proj_set, pov_set)


if __name__ == "__main__":
    main()
