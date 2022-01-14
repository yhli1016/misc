#! /usr/bin/env python
import numpy as np
import kpoints as kpt
import hamiltonian as ham
import matplotlib.pyplot as plt
import time


# Build the k-path in reciprocal space.
lattice = np.array([
[   2.458075766398899,   0.000000000000000,   0.000000000000000],
[  -1.229037883199450,   2.128755065595607,   0.000000000000000],
[   0.000000000000000,   0.000000000000000,  15.000014072326660],
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
tbmod.read_hr("graph_hr.dat", threshold=0.0)
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
plt.savefig("graph.png")
t1 = time.time()
print("Time elapsed in ploting band structure: %fs" % (t1 - t0))
