import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# Load data
fit_exp = np.array([
[ 5.900000, 147.6],
[ 7.800000, 123.5],
[10.100000,  92.1],
[14.300000,  71.4],
[26.200000,  69.5]])
fit_cal = np.loadtxt("theo_fi_1.dat")

# Global settings
plt.rc("text", usetex=True)
plt.rc("font", size=18)
fig, axes = plt.subplots(figsize=(6, 4.5))

# Plot data
axes.plot(fit_exp[:, 0], fit_exp[:, 1], "ro", label=r"\textbf{Exp.}")
axes.plot(fit_cal[:, 0], fit_cal[:, 4], "b", label=r"\textbf{Fitted}")

# Set ticks
axes.set_xlabel(r"$\mathbf{d (nm)}$")
axes.set_ylabel(r"$\mathbf{E_b (meV)}$")
axes.set_xlim([2, 30])
axes.set_ylim([40, 160])
#axes.set_xticks([i for i in range(0, 35, 5)])
axes.set_xticklabels([r"\textbf{%g}" % t for t in axes.get_xticks()])
axes.set_yticklabels([r"\textbf{%g}" % t for t in axes.get_yticks()])
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
fig.savefig("eb_fit.pdf")
