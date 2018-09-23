import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def func(k, x):
    return k[0] + k[1] * k[2] ** x


def initial_guess(x, y):
    a0 = y[-1]
    x1 = x[0]
    y1 = y[0] - a0
    x2 = x[-2]
    y2 = y[-2] - a0
    c0 = (y1 / y2) ** (1 / (x1 - x2))
    b0 = y1 / (c0 ** x1)
    k0 = np.array([a0, b0, c0])
    return k0


def residuals(k, x, y):
    return func(k, x) - y


# Load and parse data
data = np.loadtxt("fit_nseries.dat")
x = data[:, 0]
y = -data[:, 2]

# Least square fitting
k0 = initial_guess(x, y)
result = opt.leastsq(residuals, k0, args=(x, y))
print(result)

# Plot original data and fitted curve
x_plot = np.linspace(np.min(x), np.max(x), 100)
y_plot = func(result[0], x_plot)

# Global settings
plt.rc("text", usetex=True)
plt.rc("font", size=18)
fig, axes = plt.subplots(figsize=(6, 4.5))

# Plot data
axes.plot(x, y, "ro", label=r"\textbf{Calc.}")
axes.plot(x_plot, y_plot, "b", label=r"\textbf{Fitted}")

# Set ticks
axes.set_xlabel(r"$\mathbf{N_{max}}$")
axes.set_ylabel(r"$\mathbf{E_b/R}$")
axes.set_xlim([0, 220])
axes.set_ylim([14, 17])
axes.set_xticklabels([r"\textbf{%g}" % t for t in axes.get_xticks()])
axes.set_yticklabels([r"\textbf{%4.1f}" % t for t in axes.get_yticks()])
axes.minorticks_on()
axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))
axes.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
axes.tick_params(which="major", direction="in", width=1.5, length=8)
axes.tick_params(which="minor", direction="in", width=1.5, length=4)

# Set spines
for pos in ("top", "bottom", "left", "right"):
    axes.spines[pos].set_linewidth(1.5)

# Set legend
axes.legend(edgecolor="w")

# Save figure
fig.tight_layout()
fig.savefig("fit_nseries_eb.pdf")
