import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# Load and process data
dat = np.loadtxt("eps2.dat")
xco = dat[:, 0].reshape((10, 10)).T
yco = dat[:, 2].reshape((10, 10)).T
zco = dat[:, 3].reshape((10, 10)).T
eps1 = np.array([7.0 for i in range(5)])
eps2 = np.array([3.74 for i in range(5)])
d = np.array([5.9, 7.8, 10.1, 14.3, 26.2])

# Interpolate
new_func = interp.RectBivariateSpline(xco[0, :], yco[:, 0], zco)
xfi = np.linspace(xco.min(), xco.max(), 100)
yfi = np.linspace(yco.min(), yco.max(), 100)
zfi = new_func(xfi, yfi)

# Plot data
plt.rc("text", usetex=True)
plt.rc("font", size=16)
fig, axes = plt.subplots(figsize=(6, 4.5))

img = axes.contourf(xfi, yfi, zfi, 15, cmap="jet")
axes.plot(d, eps2, "wo")

# Set ticks
axes.set_title(r"\textbf{Exciton radius}")
axes.set_xlabel(r"$\mathbf{d (nm)}$")
axes.set_ylabel(r"$\mathbf{\epsilon_2}$")
axes.set_xlim([1.5, 30])
axes.set_xticks([i for i in range(5, 35, 5)])
axes.set_xticklabels([r"\textbf{%g}" % t for t in axes.get_xticks()])
axes.set_yticklabels([r"\textbf{%g}" % t for t in axes.get_yticks()])
axes.minorticks_on()
axes.xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
axes.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
axes.tick_params(which="major", direction="in", width=1.5, length=8)
axes.tick_params(which="minor", direction="in", width=1.5, length=4)

# Set spines
for pos in ("top", "bottom", "left", "right"):
    axes.spines[pos].set_linewidth(1.5)

# Set colorbar
fig.colorbar(img, label=r"\textbf{nm}", format=r"\textbf{%3.1f}")

# Save figure
fig.tight_layout()
fig.savefig("a2.pdf")
