"""
Some papers say that if there is a vacuum layer between the substrate and the
perovskite layer, than the bottom barrier of the quantum well should be also
vacuum (eps1=1). This script checks if this statement yields better exciton
binding energies and radii.
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


def residuals(x):
    global eps2, a_B, R, k, q, l, zgrid, sin_zgrid, num_series
    eb_fit_list = np.zeros(eb_ref_list.shape)

    # Update epsilon_2
    x = comm.bcast(x, root=0)
    eps2 = x
    for i, qw_width in enumerate(qw_width_list):
        # Update global variables
        a_B = eps2 / mu
        R = mu / (2 * eps2**2)
        k = (eps2 - eps3) / (eps2 + eps3)
        q = (eps2 - eps1) / (eps2 + eps1) * k
        l = qw_width / bohr2ang / 2
        zgrid = np.linspace(-2*l, 0, num_grid)
        sin_zgrid = np.sin(np.pi * zgrid / (2 * l))**2
        if abs(q) <= 1.0e-4:
            num_series = 10
        else:
            num_series = int(log(series_cutoff, abs(q)))

        # Calculate residual
        xopt, fopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=1.0, xtol=1.0e-7,
                                                 ftol=1.0e-7, full_output=True, disp=False)
        eb_meV = -fopt * R * har2eV * 1000
        eb_fit_list[i] = eb_meV
    return eb_fit_list - eb_ref_list


# Physical constants
bohr2ang = 0.52917721092
har2eV = 27.211385

# Parameters of the quantum well
#
# mu: reduced effective mass of the electron-hole pair
# qw_width: width of the quantum well in angstroms
# eb_ref: exciton binding energy from experiments in meV
# eps1: relative dielectric constant of the bottom barrier
# eps2: relative dielectric constant of the quantum well
# eps3: relative dielectric constant of the top barrier
mu = 0.117
qw_width_list = [59, 78, 101, 143, 262]
eb_ref_list = np.array([147.6, 123.5, 92.1, 71.4, 69.5])
eps1 = 1.0
eps2 = None
eps3 = 1.0

# Parameters controlling the integral
#
# num_grid: number of discrete grid points to evaluate the Coulomb integral
#           over z-axis
# num_series: number of series to expand the Coulomb integral
# series_cutoff: cutoff for determining num_series
num_grid = 100
num_series = None
series_cutoff = 1.0e-4

# Initialize MPI environment
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
jobs = [[i for i in range(num_grid) if i % size == k] for k in range(size)]

# Now we begin with the fitting.
result = opt.leastsq(residuals, x0=4.8, xtol=1.0e-3, ftol=1.0e-3)
if rank == 0:
   print(result)
