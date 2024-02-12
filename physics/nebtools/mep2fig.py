#! /usr/bin/env python
"""Plot mep to figures."""

from typing import List

from ase import Atoms
from ase.io import read, write

import nebtools.neb as neb
from nebtools.render import run_pov


class Render:
    def __init__(self, images: List[Atoms]):
        self.images = images

        # Normalization parameters
        self.correct_pbc = True
        self.align_image = True
        self.center_image = False
        self.rotate = False
        self.center_z = 0.2

        # Plotting parameters
        # Allowed styles: ase2 ase3 glass simple pale intermediate vmd jmol
        self.backend = "ase"
        self.style = "jmol"
        self.figure_width = 800
        self.figure_names = [str(_) for _ in range(len(images))]

        # Documentation of the parameters:
        # https://wiki.fysik.dtu.dk/ase/ase/io/io.html?highlight=pov
        num_atom_max = max([len(_) for _ in self.images])
        self.proj_set = {
            'show_unit_cell': 0,
        }
        self.pov_set = {
            'transparent': True,
            'background': 'White',
            'textures': [self.style for _ in range(num_atom_max)],
        }

    def run(self) -> None:
        """Common interface to render a list of images."""
        # Correct pbc and align images
        ref_image = self.images[0]
        for i, image in enumerate(self.images):
            if self.correct_pbc and i > 0:
                neb.correct_pbc(ref_image, image)
            if self.align_image and i > 0:
                neb.align_image(ref_image, image)
            if self.center_image:
                neb.center_image(image, self.center_z)

        # Render images
        for i, image in enumerate(self.images):
            if self.rotate:
                image.rotate(-90, "x", rotate_cell=True)
            if self.backend == "ase":
                write(f"{self.figure_names[i]}.png", image)
            else:
                run_pov(f"{self.figure_names[i]}", image,
                        **self.proj_set, **self.pov_set)


def main():
    """Render a mep from trajectory file."""
    images = read("mep.traj", index=":")
    render = Render(images)
    render.run()


def main2():
    """Render a list of separate images. """
    # Load images
    tags = ["co", "cooh", "hcoo", "ts1co", "ts1cooh", "ts1hcoo"]
    species = ["Cu", "Fe"]
    images, prefixes = [], []
    for t in tags:
        for s in species:
            images.append(read(f"CONTCAR_{t}-{s}", index="-1"))
            prefixes.append(f"{t}-{s}")
    render = Render(images)

    # Render the images
    # Top view
    render.rotate = False
    render.figure_names = [f"{_}.top" for _ in prefixes]
    render.run()

    # Side view
    render.rotate = True
    render.figure_names = [f"{_}.side" for _ in prefixes]
    render.run()


if __name__ == "__main__":
    main()
