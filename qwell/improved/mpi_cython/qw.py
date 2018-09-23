"""
Boosted version based on cython. Tests show that an enhancement of 76% can be achieved.

original version:
    3.050000    3.735999   69.374282
real	4m51.240s
user	9m30.697s
sys	0m0.763s

cython version:
    3.050000    3.735999   69.374282
real	1m9.137s
user	2m14.989s
sys	0m0.192s
"""
import numpy as np
import scipy.optimize as opt
from mpi4py import MPI
from core import integ_zh


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
mu = 0.117
qw_width = 262
eps1 = 8.1
eps2 = 3.05
eps3 = 1.0

# Parameters controlling the integral
#
# num_grid: number of discrete grid points to evaluate the Coulomb integral
#           over z-axis
# num_series: number of series to expand the Coulomb integral
num_grid = 100
num_series = 20

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

# Now we begin with the variation.
xopt, fopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=1.0, xtol=1.0e-7,
                                         ftol=1.0e-7, full_output=True, disp=False)
if rank == 0:
    a_nm = xopt[0] * a_B * bohr2ang / 10
    eb_meV = -fopt * R * har2eV * 1000
    print("%12.6f%12.6f%12.6f" % (eps2, a_nm, eb_meV))
