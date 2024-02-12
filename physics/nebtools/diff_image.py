#! /usr/bin/env python
"""Evaluate the differences between given image and reference mep."""

import sys

from ase.io import read

from nebtools.neb import align_image, diff_image


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
