#! /usr/bin/env python
"""Get the maximum force for each step during ionic relaxation."""

import numpy as np
from ase.io import write
from ase.io.vasp import read_vasp_out

from nebtools.idpp import norm

images = read_vasp_out("OUTCAR", index=":")

# Get positions and forces
for i, image in enumerate(images):
	forces = image.get_forces()
	max_force = max([norm(f) for f in forces])
	print("%4d%8.2f" % (i, max_force))
