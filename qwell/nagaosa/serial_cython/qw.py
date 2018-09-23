"""
A boosted version of qw_serial based on Cython. Tests show that an enhancement
of 56% can be achieved.

original version:
0.568943786621 -6.09682972195
real    17m6.970s
user    17m5.763s
sys 0m0.049s

cython version:
0.570206642151 -6.08242664268
real    7m32.140s
user    7m28.952s
sys     0m0.117s
"""
import numpy as np
import scipy.optimize as opt
from core import calc_eb


# Parameters of the quantum well
#
# mu: reduced effective mass of the electron-hole pair
# qw_width: width of the quantum well in angstroms (2*l)
# epsilon_1: relative dielectric constant of the quantum well
# epsilon_2: relative dielectric constant of the barrier
mu = 0.06826386021
qw_width = 40
epsilon_1 = 12.9
epsilon_2 = 5.6

# Parameters controlling the integral
#
# num_grid: number of discrete grid points to evaluate the Coulomb integral
#           over z-axis
# num_series: number of series to expand the Coulomb integral
num_grid = 100
num_series = 20

# Physical constants
ang2bohr = 1.88972612457

# Derived variables used throughout the code
a_B = epsilon_1 / mu
R = mu / (2 * epsilon_1**2)
q = (epsilon_1 - epsilon_2) / (epsilon_1 + epsilon_2)
l = qw_width * ang2bohr / 2
zgrid = np.linspace(-l, l, num_grid)
cos_zgrid = np.cos(np.pi * zgrid / (2 * l))**2

# Now we begin with the variation.
args = (mu, epsilon_1, epsilon_2, num_grid, num_series, a_B, R, q, l, zgrid, cos_zgrid)
xopt, yopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=1.0, args=args, xtol=1.0e-7,
                                         ftol=1.0e-7, full_output=True, disp=False)
print(xopt[0], yopt)
