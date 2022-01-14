#! /usr/bin/env python
import numpy as np
import kpoints as kpt
import hamiltonian as ham
import matplotlib.pyplot as plt
import time


# Build the k-path in reciprocal space.
lattice = np.array([
[   3.181400000,   0.000000000,   0.000000000],
[  -1.590690383,   2.755178772,   0.000000000],
[   0.000000000,   0.000000000,  15.900000000],
])
hs_kpoints = np.array([
    [0.0, 0.0, 0.0],
    [1./3, 1./3, 0.0],
    [1./2, 0.0, 0.0],
    [0.0, 0.0, 0.0],
])
kpath = kpt.gen_kpath(hs_kpoints, [40, 40, 40])
kdist = kpt.gen_kdist(lattice, kpath)

# Set up the Hamiltonian.
t0 = time.time()
tbmod = ham.TBModel()
tbmod.read_hr("mos2_hr.dat", threshold=0.0)
t1 = time.time()
print("Time elapsed in setting up Hamiltonian: %fs" % (t1 - t0))

# Calculate projection-resolved band structure.
t0 = time.time()
energies, projection = tbmod.eval_energies(kpath)
t1 = time.time()
print("Time elapsed in calculating band structure: %fs" % (t1 - t0))

# Plot the band structure.
t0 = time.time()
energies_sorted = np.array([np.sort(row) for row in energies])
for i in range(energies.shape[1]):
    plt.scatter(kdist, energies[:, i], s=projection[:, i]*10, c="r")
    plt.plot(kdist, energies_sorted[:, i], color="gray", linewidth=0.5)
plt.savefig("mos2.png")
t1 = time.time()
print("Time elapsed in ploting band structure: %fs" % (t1 - t0))
