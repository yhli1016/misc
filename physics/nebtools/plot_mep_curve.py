#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def cubic_fit(xco, yco):
    """
    We build a cubic polynomial that satisfies
        f(x0) = y0
        f(x1) = y1
        f'(x0) = 0
        f'(x1) = 0.
    """
    x0, y0 = xco[0], yco[0]
    x1, y1 = xco[1], yco[1]
    A = np.array([
            [x0**3,   x0**2, x0, 1],
            [x1**3,   x1**2, x1, 1],
            [3*x0**2,  2*x0,  1, 0],
            [3*x1**2,  2*x1,  1, 0]])
    b = np.array([y0, y1, 0, 0])
    coeff = np.linalg.solve(A, b)
    xfi = np.linspace(x0, x1, 51)
    yfi = np.polyval(coeff, xfi)
    return xfi, yfi


def add_mep(x, y, label_color, **kwargs):
    assert len(x) == len(y) == len(label_color)
    react_coord = np.array(x)
    energy = np.array(y)
    energy = energy - energy[0]
    # Draw the line connecting all images
    react_coord_fi = np.array([])
    energy_fi = np.array([])
    for i in range(energy.size-1):
        xfi, yfi = cubic_fit(react_coord[i:i+2], energy[i:i+2])
        react_coord_fi = np.append(react_coord_fi, xfi)
        energy_fi = np.append(energy_fi, yfi)
    axes.plot(react_coord_fi, energy_fi, linewidth=line_width, **kwargs)
    # Add energy levels
    for i, coord in enumerate(react_coord):
        label_x = np.array([coord - label_length / 2, coord + label_length / 2])
        label_y = np.array([energy[i], energy[i]])
        axes.plot(label_x, label_y, linewidth=label_width, color=label_color[i])


# Size and resolution of the figure in inches (1 inch = 2.54 cm)
figure_size = (5, 4)
figure_dpi = 300

# Font settings
font_size = 12
font_family = "Arial"
font_weight = "bold"

# Axes and spines settings
ymin = -2.0
ymax = 2.0
axes_width = 1.5
spines_width = 1.5

# Length, width of lines and labels
line_width = 1.0
label_length = 0.3
label_width = 1.5

# Global settings and create the figure
plt.rc("font", size=font_size, family=font_family, weight=font_weight)
fig, axes = plt.subplots(figsize=figure_size)

# Plot the MEPs
# Names for colors: (b)lue, (r)ed, (g)reen, (c)yan, (m)agenta, blac(k), (w)hite
# Allowd linestyles are "-", "--", "-.", ":"
x = [0, 1, 2, 3, 4]
add_mep(x, [0, 3.04, 2.50, 4.34, 2.67],   ["k", "b", "r", "c" ,"m"], color="k", linestyle="-", label="Path-A")
add_mep(x, [0, 2.70, 2.95, 1.35, 2.00],   ["k", "b", "r", "c" ,"m"], color="b", linestyle="-", label="Path-B")

# Adjusting the figures
axes.set_xlabel("Reaction coordinate", fontsize="large", weight=font_weight)
axes.set_ylabel("Energy (eV)", fontsize="large", weight=font_weight)

axes.set_xticks([])
#axes.set_ylim(ymin, ymax)
axes.tick_params(axis="y", width=axes_width)

# Hide or show borders
for key in ("top", "bottom", "right"):
    #axes.spines[key].set_visible(False)
    axes.spines[key].set_linewidth(spines_width)

# Left border is always shown
axes.spines["left"].set_linewidth(spines_width)

# Legend
axes.legend(edgecolor="w")

# Save the figure
fig.tight_layout()
fig.savefig("mep_curve.png", dpi=figure_dpi)
