import numpy as np
import matplotlib.pyplot as plt


prefix = "eps1_8.1_" + "25"

# Load data
fit_exp = np.loadtxt(prefix + "/fit_exp.dat")
fit_cal = np.loadtxt(prefix + "/fit_theo_fi.dat")

# Plot the figure
plt.rc("font", size=12)
fig, axes = plt.subplots()

axes.plot(fit_exp[:, 0], fit_exp[:, 1], "o", label="Exp.", color="r")
axes.plot(fit_cal[:, 0], fit_cal[:, 3], label="Theory", color="b")

axes.set_xlim([0, 30])
axes.set_ylim([50, 200])
axes.set_xlabel("Thickness (nm)", fontsize="large")
axes.set_ylabel("Binding energy (meV)", fontsize="large")
axes.tick_params(which="both", direction="in")
axes.minorticks_on()
axes.legend(edgecolor="w")
plt.savefig(prefix + "/fit.pdf")
