#! /usr/bin/env python
"""Convert CONTCARs of NEB output to traj."""

from ase.io import read, write

from check_conv import get_int


def main():
    num_image = get_int("INCAR", "IMAGES")
    images = []
    for i in range(num_image+2):
        try:
            image = read("%02d/CONTCAR"%i, index=-1)
        except FileNotFoundError:
            image = read("%02d/POSCAR"%i, index=-1)
        finally:
            images.append(image)
    write("mep.traj", images)


if __name__ == "__main__":
    main()
