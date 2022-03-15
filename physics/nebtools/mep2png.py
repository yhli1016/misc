#! /usr/bin/env python
"""Plot mep to figures in png format."""

from ase.io import read, write

import idpp


def main():
    # TODO: implement 'center_image'
    correct_pbc = True
    align_image = True
    rotate = False
    selected_images = "all"

    # Load and select images
    images = read("mep.traj", index=":")
    ref_image = images[0]
    if selected_images == "all":
        selected_images = range(len(images))
    images = [images[_] for _ in selected_images]

    # Correct pbc and align images
    for image in images[1:]:
        if correct_pbc:
            idpp.correct_pbc(ref_image, image, selected_atoms="all")
        if align_image:
            idpp.align_image(ref_image, image, selected_atoms="all")

    # Save to figure
    for i, image in enumerate(images):
        if rotate:
            image.rotate(-90, "x", rotate_cell=True)
        write(f"{selected_images[i]}.png", image)


if __name__ == "__main__":
    main()
