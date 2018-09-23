"""
Demonstrate the convergence of exciton binding energy on num_series, using data
point E in nagaosa's paper.
"""
import numpy as np
import scipy.optimize as opt
from mpi4py import MPI
from core import integ_zh


def calc_eb(x):
    # Convert the unit of x from exciton Bohr radius to absolute Bohr
    a = x * a_B
    a = comm.bcast(a, root=0)

    # Evaluation the kinetic term
    kinetic_energy = 1.0 / (2.0 * mu * a**2)

    # Evaluate the Coulomb attraction term
    integral_ze = 0.0
    for i in jobs[rank]:
        ze = zgrid.item(i)
        integral_zh = integ_zh(ze, a, num_series, q, l, zgrid, cos_zgrid)
        integral_ze += cos_zgrid.item(i) * integral_zh
    integral_ze = comm.allreduce(integral_ze)
    dS = (zgrid.item(1) - zgrid.item(0)) ** 2
    potential_energy = -2.0 / (epsilon_1 * l**2 * a) * integral_ze * dS

    # Evaluate the binding energy in effective Rydberg energy
    binding_energy = (kinetic_energy + potential_energy) / R
    return binding_energy


# Physical constants
ang2bohr = 1.88972612457

# Parameters of the quantum well
#
# mu: reduced effective mass of the electron-hole pair
# qw_width: width of the quantum well in angstroms (2*l)
# epsilon_1: relative dielectric constant of the quantum well
# epsilon_2: relative dielectric constant of the barrier
mu = 0.06826386021
qw_width = 100
epsilon_1 = 12.9
epsilon_2 = epsilon_1 / 100

# Parameters controlling the integral
#
# num_grid: number of discrete grid points to evaluate the Coulomb integral
#           over z-axis
# num_series: number of series to expand the Coulomb integral
num_grid = 100
num_series = None

# Derived variables used throughout the code
a_B = epsilon_1 / mu
R = mu / (2 * epsilon_1**2)
q = (epsilon_1 - epsilon_2) / (epsilon_1 + epsilon_2)
l = qw_width * ang2bohr / 2
zgrid = np.linspace(-l, l, num_grid)
cos_zgrid = np.cos(np.pi * zgrid / (2 * l))**2

# Initialize MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
jobs = [[i for i in range(num_grid) if i % size == k] for k in range(size)]

# Main loop
x0 = 1.0
for num_series in range(20, 220, 20):
    xopt, yopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=x0,
                                                      xtol=1.0e-7,
                                                      ftol=1.0e-7,
                                                      full_output=True,
                                                      disp=False)
    if rank == 0:
        print(num_series, xopt[0], yopt, flush=True)
    x0 = xopt[0]
