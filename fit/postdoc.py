#! /usr/bin/env python
import numpy as np
from scipy import optimize as opt


def obj_func(k):
    E0 = k[0]
    alpha = k[1]
    beta = k[2]
    gamma = k[3]
    delta = k[4]
    Npi0 = k[5]
    return E0 + alpha * N2c / NT + beta * N3c / NT + gamma * N4c / NT \
              + delta * ((Npi - Npi0) / NT)**2


def residual(k):
    return obj_func(k) - Ec


# Load data
dat = np.loadtxt("EC_this.dat")
N2c = dat[:, 0]
N3c = dat[:, 1]
N4c = dat[:, 2]
Npi = dat[:, 3]
NT = dat[:, 4]
Ec = dat[:, 5]

# Least-square fit
k0 = [5.0, 0.5, -0.1, 0.3, 0.5, 0.2]
result = opt.leastsq(residual, k0)
print("optimal parameters:", result[0])
print("error code:", result[1])

# Output to file
N2c = 6
N3c = None
N4c = 8
Npi = None
NT = 36

k_opt = result[0]
x_fine = np.linspace(0, NT - N2c - N4c, 100)
for x in x_fine:
    N3c = x
    Npi = NT - N2c - N4c - N3c
    print("%12.5f%12.5f" % (N3c, obj_func(k_opt)))
