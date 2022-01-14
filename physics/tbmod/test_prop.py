#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import hamiltonian as ham


tbmod = ham.TBModel()
tbmod.read_hr("hr.dat")

eng1, vec = tbmod.eval_energies(np.array([[0.5, 0.0, 0.0]]))
eng2, dos, peaks = tbmod.propagate(np.array([0.5, 0.0, 0.0]), 0.5, 1000000)

plt.plot(eng2, dos)
plt.show()
print(eng1)
print(peaks)
