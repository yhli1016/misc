import numpy as np
import kpoints as kpt
import hamiltonian as ham
import matplotlib.pyplot as plt

lattice = np.array([
    [2.027733378, 3.512137236, 0.000000000],
    [-2.027733378, 3.512137236, 0.000000000],
    [0.000000000, 0.000000000, 20.500000000],
])

hs_kpoints = np.array([
    [0.0, 0.0, 0.0],
    [2/3, 1/3, 0.0],
    [1/2, 0.0, 0.0],
    [0.0, 0.0, 0.0],
])
kpath = kpt.gen_kpath(hs_kpoints, [40, 40, 40])
kdist = kpt.gen_kdist(lattice, kpath)

tbmod = ham.TBModel(16)
tbmod.read_hr("hr.dat")

energies, projection = tbmod.eval_energies(kpath, [9, 13])
energies_sorted = np.array([np.sort(row) for row in energies])
for i in range(energies.shape[1]):
    plt.scatter(kdist, energies[:, i], s=projection[:, i] * 10, c="r")
    plt.plot(kdist, energies_sorted[:, i], color="gray", linewidth=0.5)
plt.show()
