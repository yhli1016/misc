#! /usr/bin/env python
"""Plot mep to figures in png format."""

from ase.io import read, write

import idpp


def main():
    # TODO: implement 'center_image'
    correct_pbc = True
    align_image = True
    rotate = True

    images = read("mep.traj", index=":")

    # Correct pbc and align images
    ref_image = images[0]
    for image in images[1:]:
        if correct_pbc:
            idpp.correct_pbc(ref_image, image, selected_atoms="all")
        if align_image:
            idpp.align_image(ref_image, image, selected_atoms="all")

    # Save to figure
    for i, image in enumerate(images):
        if rotate:
            image.rotate(-90, "x", rotate_cell=True)
        write(f"{i}.png", image)


if __name__ == "__main__":
    main()
