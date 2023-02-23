#! /usr/bin/env python
"""Plot mep to figures."""

import os
import sys
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
    figure_height = int(figure_width * aspect_ratio)
    config.set("default", "Width", str(figure_width))
    config.set("default", "Height", str(figure_height))

    # Update ini file
    with open(f"{prefix}.ini", "w") as ini_file:
        config.write(ini_file)


def run_pov(prefix, image, figure_width=800, proj_set=None, pov_set=None):
    """
    Run POV-Ray to render an image.

    :param str prefix: prefix prepended to output files
    :param ase.Atoms image: image to render
    :param int figure_width: width of figure in pixel
    :param dict proj_set: general projection settings
    :param dict pov_set: POV-Ray settings
    :return: None
    """
    # Write ini and pov file
    write(f"{prefix}.pov", image, **proj_set, povray_settings=pov_set)
    adjust_pov_ini(prefix, figure_width)

    # Run POV-Ray
    if sys.platform == "linux":
        os.system(f"povray {prefix}[default]")
    else:
        os.system(f"pvengine64 /exit {prefix}[default]")

    # Clean up
    for suffix in ("ini", "pov"):
        os.remove(f"{prefix}.{suffix}")


def render_images(images, correct_pbc=False, align_image=False,
                  center_image=False, center_z=0.5, rotate=False,
                  render="ase", figure_names=None, **kwargs):
    """Common interface to render a list of images."""
    # Correct pbc and align images
    ref_image = images[0]
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
        if render == "ase":
            write(f"{figure_names[i]}.png", image)
        else:
            run_pov(f"{figure_names[i]}", image, **kwargs)


def main():
    """Render a mep from trajectory file."""
    # Normalization parameters
    correct_pbc = True
    align_image = True
    center_image = False
    rotate = False
    center_z = 0.2

    # Plotting parameters
    # Allowed styles: ase2 ase3 glass simple pale intermediate vmd jmol
    selected_images = "all"
    render = "ase"
    style = "jmol"
    figure_width = 800

    # Load images
    images = read("mep.traj", index=":")
    if selected_images == "all":
        selected_images = range(len(images))
    images = [images[_] for _ in selected_images]

    # Documentation of the parameters:
    # https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=pov
    num_atom = [len(_) for _ in images]
    textures = [style for _ in range(max(num_atom))]
    proj_set = {
        'show_unit_cell': 0,
    }
    pov_set = {
        'transparent': True,
        'background': 'White',
        'textures': textures,
    }

    # Render the images
    render_images(images, correct_pbc=correct_pbc, align_image=align_image,
                  center_image=center_image, center_z=center_z, rotate=rotate,
                  render=render, figure_width=figure_width,
                  figure_names=selected_images, proj_set=proj_set,
                  pov_set=pov_set)


def main2():
    """Render a list of separate images. """
    # Normalization parameters
    correct_pbc = False
    align_image = False
    center_image = False
    center_z = 0.2

    # Plotting parameters
    # Allowed styles: ase2 ase3 glass simple pale intermediate vmd jmol
    render = "ase"
    style = "jmol"
    figure_width = 800

    # Load images
    tags = ["co", "cooh", "hcoo", "ts1co", "ts1cooh", "ts1hcoo"]
    species = ["Cu", "Fe"]
    images, prefixes = [], []
    for t in tags:
        for s in species:
            images.append(read(f"CONTCAR_{t}-{s}", index="-1"))
            prefixes.append(f"{t}-{s}")

    # Documentation of the parameters:
    # https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=pov
    num_atom = [len(_) for _ in images]
    textures = [style for _ in range(max(num_atom))]
    proj_set = {
        'show_unit_cell': 0,
    }
    pov_set = {
        'transparent': True,
        'background': 'White',
        'textures': textures,
    }

    # Render the images
    # Top view
    rotate = False
    figure_names = [f"{_}.top" for _ in prefixes]
    render_images(images, correct_pbc=correct_pbc, align_image=align_image,
                  center_image=center_image, center_z=center_z, rotate=rotate,
                  render=render, figure_width=figure_width,
                  figure_names=figure_names, proj_set=proj_set,
                  pov_set=pov_set)
    # Side view
    rotate = True
    figure_names = [f"{_}.side" for _ in prefixes]
    render_images(images, correct_pbc=correct_pbc, align_image=align_image,
                  center_image=center_image, center_z=center_z, rotate=rotate,
                  render=render, figure_width=figure_width,
                  figure_names=figure_names, proj_set=proj_set,
                  pov_set=pov_set)


if __name__ == "__main__":
    main()
