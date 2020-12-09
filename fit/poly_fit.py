#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


# Load and parse data
data = np.loadtxt("fit.dat")
x = data[:, 0]
y = data[:, 1]

# Least square fitting
p = np.poly1d(np.polyfit(x, y, 5))
q = np.roots(np.polyder(p))
print(q)

# Plot original data and fitted curve
x_plot = np.linspace(np.min(x), np.max(x), 100)
y_plot = p(x_plot)
plt.plot(x, y, "o", color="r")
plt.plot(x_plot, y_plot, color="b")
plt.savefig("fit.png")
