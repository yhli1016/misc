#! /usr/bin/env python
import numpy as np
import kpoints as kpt
import hamiltonian as ham
import matplotlib.pyplot as plt
import time


# Build the k-path in reciprocal space.
lattice = np.array([
    [2.027733378, 3.512137236, 0.000000000],
    [-2.027733378, 3.512137236, 0.000000000],
    [0.000000000, 0.000000000, 20.500000000],
])
hs_kpoints = np.array([
    [0.0, 0.0, 0.0],
    [2./3, 1./3, 0.0],
    [1./2, 0.0, 0.0],
    [0.0, 0.0, 0.0],
])
kpath = kpt.gen_kpath(hs_kpoints, [40, 40, 40])
kdist = kpt.gen_kdist(lattice, kpath)

# Set up the Hamiltonian.
t0 = time.time()
tbmod = ham.TBModel()
tbmod.read_hr("wannier/inse_hr.dat")
t1 = time.time()
print("Time elapsed in setting up Hamiltonian: %fs" % (t1 - t0))

# Calculate band structure by time propagation.
t0 = time.time()
dt = 1.0
nstep = 10000
peaks_total = []
for ik, kpt in enumerate(kpath):
    print("Dealing with kpt #%s" % ik)
    eng, dos, peaks = tbmod.propagate(kpt, dt, nstep)
    peaks_total.append(peaks)
t1 = time.time()
print("Time elapsed in calculating band structure: %fs" % (t1 - t0))

# Plot the band structure.
t0 = time.time()
for ik, klen in enumerate(kdist):
    peaks = peaks_total[ik]
    x = [klen for _ in range(peaks.shape[0])]
    plt.scatter(x, peaks, s=2, c='r')
plt.show()
t1 = time.time()
print("Time elapsed in ploting band structure: %fs" % (t1 - t0))
