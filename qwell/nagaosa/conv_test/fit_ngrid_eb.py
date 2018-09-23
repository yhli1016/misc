#! /usr/bin/env python

import numpy as np
from scipy import optimize as opt
import matplotlib.pyplot as plt


def func(k, x):
    return k[0] + k[1] / (k[2] + x)


def initial_guess(x, y):
    a0 = y[-1]
    x1 = x[0]
    y1 = y[0] - a0
    x2 = x[-2]
    y2 = y[-2] - a0
    c0 = (x1 * y1 - x2 * y2) / (y2 - y1)
    b0 = y1 * (c0 + x1)
    k0 = np.array([a0, b0, c0])
    return k0


def residuals(k, x, y):
    return func(k, x) - y


# Load and parse data
data = np.loadtxt("fit_ngrid.dat")
x = data[:, 0]
y = data[:, 4]

# Least square fitting
k0 = initial_guess(x, y)
result = opt.leastsq(residuals, k0, args=(x, y))
print("optimal parameters:", result[0])
print("error code:", result[1])

# Plot original data and fitted curve
x_plot = np.linspace(np.min(x), np.max(x), 100)
y_plot = func(result[0], x_plot)
plt.plot(x, y, "o", color="r")
plt.plot(x_plot, y_plot, color="b")
plt.savefig("fit_ngrid_eb.png")
