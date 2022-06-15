#! /usr/bin/env python
"""Plot mep to figures."""

import os
from configparser import ConfigParser

from ase.io import read, write

from nebtools import idpp


def adjust_pov_ini(prefix, figure_width=800):
    """
    Adjust POV-Ray settings in ini file.

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


def render_pov(prefix, image, image_width, proj_set, pov_set):
    """
    Render atoms via POV-Ray.

    :param str prefix: prefix prepended to output files
    :param ase.Atoms image: image to render
    :param int image_width: width of image in pixel
    :param dict proj_set: general projection settings
    :param dict pov_set: POV-Ray settings
    :return: None
    """
    # Write ini and pov file
    write(f"{prefix}.pov", image, **proj_set, povray_settings=pov_set)
    adjust_pov_ini(prefix, image_width)

    # Run POV-Ray
    # os.system(f"pvengine64 /exit {prefix}[default]")
    os.system(f"povray {prefix}[default]")

    # Clean up
    for suffix in ("ini", "pov"):
        os.remove(f"{prefix}.{suffix}")


def main():
    # Normalization parameters
    correct_pbc = True
    align_image = True
    center_image = False
    rotate = False
    center_z = 0.2

    # Plotting parameters
    selected_images = "all"
    fig_format = "png"
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
    for i, image in enumerate(images):
        if correct_pbc and i > 0:
            idpp.correct_pbc(ref_image, image)
        if align_image and i > 0:
            idpp.align_image(ref_image, image)
        if center_image:
            idpp.center_image(image, center_z)

    # Render images
    for i, image in enumerate(images):
        if rotate:
            image.rotate(-90, "x", rotate_cell=True)
        if fig_format == "png":
            write(f"{selected_images[i]}.png", image)
        else:
            render_pov(f"{selected_images[i]}", image, figure_width,
                       proj_set, pov_set)


if __name__ == "__main__":
    main()
