import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


prefix = "eps1_8.1" + "_10"

# Load and process data
dat = np.loadtxt(prefix + "/eps2.dat")
xco = dat[:, 0].reshape((10, 10)).T
yco = dat[:, 2].reshape((10, 10)).T
zco = dat[:, 3].reshape((10, 10)).T

# Interpolate
new_func = interp.RectBivariateSpline(xco[0, :], yco[:, 0], zco)
xfi = np.linspace(xco.min(), xco.max(), 100)
yfi = np.linspace(yco.min(), yco.max(), 100)
zfi = new_func(xfi, yfi)

# Plot
plt.rc("font", size=12, family="Times New Roman", weight="bold")
fig, axes = plt.subplots()

img = axes.contourf(xfi, yfi, zfi, 15, cmap="jet")

axes.set_title("Exciton radius", fontsize="large", weight="bold")
axes.set_xlabel("Thickness (nm)", fontsize="large", weight="bold")
axes.set_ylabel("$\epsilon_2$", fontsize="large", weight="bold")
axes.minorticks_on()
axes.xaxis.set_minor_locator(ticker.MultipleLocator(2.5))
axes.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
axes.tick_params(which="both", direction="in")
fig.colorbar(img, label="nm")

# Add the samples
eps1 = np.array([8.1, 8.1, 8.1, 8.1, 8.1])
eps2 = np.array([3.52, 3.52, 3.52, 3.52, 3.52])
d = np.array([5.8, 7.8, 10.1, 14.3, 26.2])
axes.plot(d, eps2, "wo")

# Save to file
fig.savefig(prefix + "/radius2.pdf")
with open (prefix + "/a_eps2.dat", "w") as f:
    for i, x in enumerate(xfi):
        for j, y in enumerate(yfi):
            f.write("%12.6f%12.6f%12.6f\n" % (x, y, zfi[j,i]))
