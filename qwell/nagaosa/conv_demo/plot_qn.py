import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def func(x, ratio):
    return (1 - 2 / (ratio + 1))**x


x = np.linspace(1, 200, 100)
y1 = func(x, 10)
y2 = func(x, 100)


# Global settings
plt.rc("text", usetex=True)
plt.rc("font", size=18)
fig, axes = plt.subplots(figsize=(6, 4.5))

# Plot data
axes.plot(x, y1, "r", label="$\mathbf{\epsilon_2 / \epsilon_1=10}$")
axes.plot(x, y2, "b", label="$\mathbf{\epsilon_2 / \epsilon_1=100}$")

# Set tickes
axes.set_xlabel(r"$\mathbf{n}$")
axes.set_ylabel(r"$\mathbf{q^n}$")
axes.set_xlim([-10, 210])
axes.set_ylim([-0.05, 1.0])
axes.set_xticklabels([r"\textbf{%g}" % t for t in axes.get_xticks()])
axes.set_yticklabels([r"\textbf{%3.1f}" % t for t in axes.get_yticks()])
axes.minorticks_on()
axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))
axes.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
axes.tick_params(which="major", direction="in", width=1.5, length=8)
axes.tick_params(which="minor", direction="in", width=1.5, length=4)

# Set spines
for pos in ("top", "bottom", "left", "right"):
    axes.spines[pos].set_linewidth(1.5)

# Set legend
axes.legend(edgecolor="w", fontsize="large")

# Save figure
fig.tight_layout()
fig.savefig("plot_qn.pdf")
