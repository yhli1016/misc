#! /usr/bin/env python
"""Get the maximum force for each step during ionic relaxation."""

from ase.io.vasp import read_vasp_out

from nebtools.neb import norm

images = read_vasp_out("OUTCAR", index=":")

# Get positions and forces
for i, image in enumerate(images):
    forces = image.get_forces()
    max_force = max([norm(f) for f in forces])
    print(f"{i:4d}{max_force:8.2f}")
