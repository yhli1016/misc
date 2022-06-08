#! /usr/bin/env python
"""Adjust image according to atomic forces."""

import numpy as np
from ase.io import write
from ase.io.vasp import read_vasp_out


image = read_vasp_out("OUTCAR", index=-1)

# Get positions and forces
positions = image.get_positions()
forces = image.get_forces()

# Adjust positions according to forces
r_max = float(input("Give max. displacement in Angstrom: "))
f_max = np.max([np.linalg.norm(f) for f in forces])
for i in range(len(image)):
    dr = forces[i] / f_max * r_max
    positions[i] += dr
image.set_positions(positions)
write("CONTCAR", image, vasp5=True, direct=True)

# Write to xyz file
with open("test.xyz", "w") as xyz:
    xyz.write("%4d\n" % len(image))
    xyz.write("Generated with ASE\n")
    for i, atom in enumerate(image):
        xyz.write("%4s%14.9f%14.9f%14.9f%14.9f%14.9f%14.9f\n" % (atom.symbol,
                   positions[i, 0], positions[i, 1], positions[i, 2],
                   forces[i, 0], forces[i, 1], forces[i, 2]))
