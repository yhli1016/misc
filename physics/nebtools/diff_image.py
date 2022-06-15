#! /usr/bin/env python
"""Evaluate the differences between given image and reference mep."""

import sys

from ase.io import read

from idpp import norm, align_image


def diff_image(ref_image, chk_image):
    """
    Evaluate the difference between given image and reference.

    :param ase.Atoms ref_image: reference image
    :param ase.Atoms chk_image: image to diff
    :return: float, total difference
    """
    ref_pos = ref_image.get_positions()
    chk_pos = chk_image.get_positions()
    diff_pos = ref_pos - chk_pos
    return sum([norm(diff_pos[i_a])**2 for i_a in range(len(ref_image))])


def main():
    chk_image = read(sys.argv[1], index="-1")
    mep = read(sys.argv[2], index=":")
    _align_image = False

    differences = []
    for ref_image in mep:
        if _align_image:
            align_image(ref_image, chk_image)
        differences.append(diff_image(ref_image, chk_image))

    for diff in differences:
        print(diff)


if __name__ == "__main__":
    main()
