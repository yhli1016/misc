#! /usr/bin/env python
"""Prepare directories and POSCARs from mep.traj for initial NEB run."""

from ase.io import read

from idpp import mep2pos


def main():
    mep_name = input("Input the name of mep to resume: ")
    mep = read(mep_name, index=":")
    mep2pos(mep)


if __name__ == "__main__":
    main()
