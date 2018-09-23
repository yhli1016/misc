"""
Test the convergence of exciton binding energy on num_grid and num_series,
using data point E in nagaosa's paper.
"""
import numpy as np
import scipy.optimize as opt
from mpi4py import MPI
from core import integ_zh
from math import log


def calc_eb(x):
    # Convert the unit of x from exciton Bohr radius to absolute Bohr
    x = comm.bcast(x, root=0)
    a = x * a_B

    # Evaluation the kinetic term
    kinetic_energy = 1.0 / (2.0 * mu * a**2)

    # Evaluate the Coulomb attraction term
    integral_ze = 0.0
    for i in jobs[rank]:
        ze = zgrid.item(i)
        integral_zh = integ_zh(ze, a, num_series, q, k, l, zgrid, sin_zgrid)
        integral_ze += sin_zgrid.item(i) * integral_zh
    integral_ze = comm.allreduce(integral_ze)
    dS = (zgrid.item(1) - zgrid.item(0)) ** 2
    potential_energy = -2.0 / (eps2 * l**2 * a) * integral_ze * dS

    # Evaluate the binding energy in effective Rydberg energy
    binding_energy = (kinetic_energy + potential_energy) / R
    return binding_energy


# Physical constants
bohr2ang = 0.52917721092
har2eV = 27.211385

# Parameters of the quantum well
#
# mu: reduced effective mass of the electron-hole pair
# qw_width: width of the quantum well in angstroms (2*l)
# eps1: relative dielectric constant of the bottom barrier
# eps2: relative dielectric constant of the quantum well
# eps3: relative dielectric constant of the top barrier
mu = 0.06826386021
qw_width = 100
eps2 = 12.9
eps1 = eps3 = eps2 / 100

# Parameters controlling the integral
#
# num_grid: number of discrete grid points to evaluate the Coulomb integral
#           over z-axis
# num_series: number of series to expand the Coulomb integral
x0 = 1.0
for num_grid in range(50, 450, 50):
    # Derived variables used throughout the code
    a_B = eps2 / mu
    R = mu / (2 * eps2**2)
    k = (eps2 - eps3) / (eps2 + eps3)
    q = (eps2 - eps1) / (eps2 + eps1) * k
    l = qw_width / bohr2ang / 2
    zgrid = np.linspace(-2*l, 0, num_grid)
    sin_zgrid = np.sin(np.pi * zgrid / (2 * l))**2

    # Initialize MPI environment
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    jobs = [[i for i in range(num_grid) if i % size == k] for k in range(size)]

    for n in range(2, 10, 2):
        conv_thr = 10**(-n)
        num_series = int(log(conv_thr, q))
        xopt, fopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=x0,
                                                          xtol=1.0e-7,
                                                          ftol=1.0e-7,
                                                          full_output=True,
                                                          disp=False)
        if rank == 0:
            print(num_grid, conv_thr, num_series, xopt[0], fopt, flush=True)
        x0 = xopt[0]
