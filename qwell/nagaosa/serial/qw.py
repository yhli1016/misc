"""
NOTES
-----
[1] Hartree atomic units are used throughout the code unless otherwise stated.

[2] The Hamiltonian of exciton should have six terms according to Nagaosa's and
PRB papers. However, only the kinetic term (Tr) corresponding to the relative
motion of electron-hole pair and the Coulomb attraction term (Veh) are dependent
on the variation parameter a. Therefore we neglect all the other terms and take
the variation of exciton binding energy (Tr + Veh) instead of the total energy.
"""
import numpy as np
import scipy.special as sp
import scipy.optimize as opt


def G(a, gamma):
    if abs(gamma) < 1.0e-9:
        return 1.0
    else:
        x = 2 * abs(gamma) / a
        return x * (np.pi / 2 * (sp.struve(1, x) - sp.y1(x)) - 1)
        

def calc_eb(x):
    # Convert the unit of x from exciton Bohr radius to absolute Bohr
    a = x * a_B

    # Evaluation the kinetic term
    kinetic_energy = 1.0 / (2.0 * mu * a**2)

    # Evaluate the Coulomb attraction term
    integral_ze = 0.0
    for i, ze in enumerate(zgrid):
        integral_zh = 0.0
        for j, zh in enumerate(zgrid):
            kernel = 0.0
            for n in range(-num_series, num_series+1):
                qn = q**abs(n)
                zhn = zh * (-1)**n + 2 * n * l
                kernel += qn * G(a, ze -zhn)
            integral_zh += cos_zgrid[j] * kernel
        integral_ze += cos_zgrid[i] * integral_zh
    dS = (zgrid[1] - zgrid[0]) ** 2
    potential_energy = -2.0 / (epsilon_1 * l**2 * a) * integral_ze * dS

    # Evaluate the binding energy in effective Rydberg energy
    binding_energy = (kinetic_energy + potential_energy) / R

    return binding_energy


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
xopt, yopt, ninter, funcalls, warnflag = opt.fmin(calc_eb, x0=1.0, xtol=1.0e-7,
                                         ftol=1.0e-7, full_output=True, disp=False)
print(xopt[0], yopt)
