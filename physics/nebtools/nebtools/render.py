"""Functions and classes for rendering images using ase or pov-ray"""

import os
import sys
from configparser import ConfigParser

from ase import Atoms
from ase.io import write


def adjust_pov_ini(prefix: str, figure_width: int = 800) -> None:
    """
    Adjust POV-Ray settings in ini file.

    :param prefix: prefix of ini file
    :param figure_width: width of figure
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


def run_pov(prefix: str,
            image: Atoms,
            figure_width: int = 800,
            proj_set: dict = None,
            pov_set: dict = None) -> None:
    """
    Run POV-Ray to render an image.

    :param prefix: prefix prepended to output files
    :param ase.Atoms image: image to render
    :param figure_width: width of figure in pixel
    :param proj_set: general projection settings
    :param pov_set: POV-Ray settings
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
