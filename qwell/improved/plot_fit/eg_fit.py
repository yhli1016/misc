import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def func(k, x):
    return k[0] + k[1] * x


def residuals(k, x, y):
    return func(k, x) - y


# Prepare data
data = np.array([
    [5.916, 2.432],
    [7.823, 2.398],
    [10.316, 2.384],
    [14.274, 2.366],
    [19.846, 2.357],
    [26.151, 2.348]])

x = 1 / data[:, 0]**2
y = data[:, 1]

# Least square fitting
k0 = [1, 0]
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
axes.plot(x*1000, y, "ro", label=r"\textbf{Exp.}")
axes.plot(x_plot*1000, y_plot, "b", label=r"\textbf{Fitted}")

# Set ticks
axes.set_xlabel(r"$\mathbf{d^{-2} (10^{-3}nm^{-2})}$")
axes.set_ylabel(r"$\mathbf{E_g (eV)}$")
axes.set_xlim([0, 30])
axes.set_ylim([2.34, 2.44])
axes.set_xticks([i for i in range(0, 35, 5)])
axes.set_xticklabels([r"\textbf{%g}" % t for t in axes.get_xticks()])
axes.set_yticklabels([r"\textbf{%4.2f}" % t for t in axes.get_yticks()])
axes.minorticks_on()
#axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))
#axes.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
axes.tick_params(which="major", direction="in", width=1.5, length=8)
axes.tick_params(which="minor", direction="in", width=1.5, length=4)

# Set spines
for pos in ("top", "bottom", "left", "right"):
    axes.spines[pos].set_linewidth(1.5)

# Set legend
axes.legend(edgecolor="w")

# Save figure
fig.tight_layout()
fig.savefig("eg_fit.pdf")
