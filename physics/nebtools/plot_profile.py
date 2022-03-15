#! /usr/bin/env python
"""Plot energy profiles of reaction paths."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


# Physical constants
EV2KJMOL = 96.4916


class Profile:
    """
    Class for representing the energy profile of a reaction path.

    Attributes
    ---------
    axes: instance of 'matplotlib.pyplot.Axes' class
        axes on which the profiles will be plotted
    default_label_color: string
        default color of labels
    default_connector: integer
        order of polynomial for connecting labels (energy levels)
        in the profile
    label_length: float
        length of labels of energy levels
    label_width: float
        width of labels of energy levels
    line_width: float
        width of lines connecting labels
    react_coord: List[int]
        x-coordinates of labels
    energy: List[float]
        y-coordinates of labels (energy levels)
    label_color: List[str]
        colors for each label
    connector: list[int]
        order of polynomial for connecting each label and the next label
    """
    def __init__(self, axes, default_label_color="k", default_connector=1,
                 label_length=0.6, label_width=1.8, line_width=0.75) -> None:
        self.axes = axes
        self.default_label_color = default_label_color
        self.default_connector = default_connector
        self.label_length = label_length
        self.label_width = label_width
        self.line_width = line_width
        self.react_coord = []
        self.energy = []
        self.label_color = []
        self.connector = []

    def add_eng(self, react_coord, energy, label_color=None, connector=None):
        """
        Add an absolute energy level to the profile.

        :param int react_coord: x-coordinate of the energy level
        :param float energy: y-coordinate of the energy level
        :param str label_color: color of the energy level
        :param int connector: order of polynomial for the label
        :return: None
        """
        self.react_coord.append(react_coord)
        self.energy.append(energy)
        if label_color is None:
            label_color = self.default_label_color
        self.label_color.append(label_color)
        if connector is None:
            connector = self.default_connector
        self.connector.append(connector)

    def add_de(self, de, dx=1, **kwargs):
        """
        Add a relative energy level to the profile.

        :param float de: relative energy w.r.t. previous energy level
        :param int dx: incremental of reaction coordinate of the energy level
        :param dict kwargs: arguments to be passed to the 'add_eng' method
        :return: None
        """
        react_coord = self.curr_x + dx
        if len(self.energy) == 0:
            energy = de
        else:
            energy = self.energy[-1] + de
        self.add_eng(react_coord, energy, **kwargs)

    @property
    def curr_x(self):
        """Get the reaction coordinate of last energy level."""
        if len(self.react_coord) == 0:
            return 0.0
        else:
            return self.react_coord[-1]

    def linear_fit(self, xco, yco):
        """Linear connector."""
        x0 = xco[0] + 0.5 * self.label_length
        x1 = xco[1] - 0.5 * self.label_length
        y0, y1 = yco[0], yco[1]
        xfi = np.linspace(x0, x1, 51)
        yfi = np.linspace(y0, y1, 51)
        return xfi, yfi

    @staticmethod
    def cubic_fit(xco, yco):
        """Cubic connector."""
        x0, y0 = xco[0], yco[0]
        x1, y1 = xco[1], yco[1]
        a = np.array([
                [x0**3,   x0**2, x0, 1],
                [x1**3,   x1**2, x1, 1],
                [3*x0**2,  2*x0,  1, 0],
                [3*x1**2,  2*x1,  1, 0]])
        b = np.array([y0, y1, 0, 0])
        c = np.linalg.solve(a, b)
        xfi = np.linspace(x0, x1, 51)
        yfi = np.polyval(c, xfi)
        return xfi, yfi

    def plot(self, unit="kjm", **kwargs):
        """
        Plot energy profile onto axes.

        :param str unit: unit for energy levels
        :param dict kwargs: arguments to be passed to Axes.plot()
        :return: None
        """
        # Prepare data
        react_coord = np.array(self.react_coord)
        if unit == "kjm":
            energy = np.array(self.energy) * EV2KJMOL
        else:
            energy = np.array(self.energy)

        # Draw the connectors
        react_coord_fi = np.array([])
        energy_fi = np.array([])
        for i in range(energy.size-1):
            if self.connector[i] == 1:
                fit_func = self.linear_fit
            else:
                fit_func = self.cubic_fit
            xfi, yfi = fit_func(react_coord[i:i+2], energy[i:i+2])
            react_coord_fi = np.append(react_coord_fi, xfi)
            energy_fi = np.append(energy_fi, yfi)
        self.axes.plot(react_coord_fi, energy_fi, linewidth=self.line_width,
                       **kwargs)

        # Add energy levels
        levels = []
        for i, coord in enumerate(react_coord):
            p0 = (coord - 0.5 * self.label_length, energy[i])
            p1 = (coord + 0.5 * self.label_length, energy[i])
            levels.append((p0, p1))
        levels = LineCollection(levels, colors=self.label_color,
                                linewidth=self.label_width)
        self.axes.add_collection(levels)

    def print(self, unit="kjm"):
        """
        Print absolute energy levels of the profile.

        :param str unit: unit for the energy levels
        :return: None
        """
        if unit == "kjm":
            for data in zip(self.react_coord, self.energy):
                print("%2d%8.2f" % (data[0], data[1]*EV2KJMOL))
        else:
            for data in zip(self.react_coord, self.energy):
                print("%2d%8.2f" % (data[0], data[1]))


def main():
    # Figure Settings
    figure_size = (7, 5)
    figure_dpi = 300

    # Font settings
    font_size = 12
    # font_family = "Arial"
    font_family = "Liberation Sans"
    font_weight = "normal"

    # Axes and spines settings
    # ymin = -2.0
    # ymax = 2.0
    axes_width = 1.5
    spines_width = 1.5

    # Change global settings and create the figure
    plt.rc("font", size=font_size, family=font_family, weight=font_weight)
    fig, axes = plt.subplots(figsize=figure_size)

    # Names for colors: (b)lue, (r)ed, (g)reen, (c)yan, (m)agenta, blac(k), (w)hite
    # Allowed line styles are "-", "--", "-.", ":"

    # for 0 coverage
    path = Profile(axes, default_label_color="k")
    path.add_de(0.0)    # co2 + sub
    path.add_de(-1.26)  # co2_fe_bend
    path.add_de(0.39)   # co_form_ts
    path.add_de(-1.23)  # co_fe_o_sub
    path.add_de(1.56)   # co + o_sub
    path.add_de(-0.24)  # h2_fe
    path.add_de(0.41)   # h2o_from1_ts
    path.add_de(-0.31)  # h_fe_h_o
    path.add_de(1.18)   # h2o_form2_ts
    path.add_de(-0.49)  # h2o_sub
    path.add_de(0.69)   # h2o + sub
    plot_args = {"color": "k", "linestyle": "--", "label": "0 coverage"}
    path.plot(**plot_args)
    path.print()

    # for 0.3 coverage
    path = Profile(axes, default_label_color="r")
    path.add_de(0.0)    # co2 + sub
    path.add_de(-1.40)  # co2_fe_bend
    path.add_de(0.92)   # co_form_ts
    path.add_de(-1.57)  # co_fe_o_sub
    path.add_de(1.76)   # co + o_sub
    path.add_de(-0.59)  # h2_fe
    path.add_de(0.56)   # h2o_from1_ts
    path.add_de(-0.92)  # h_fe_h_o
    # path.add_dE()  # h2o_form2_ts
    path.add_de(dx=2, de=1.13)  # h2o_sub
    path.add_de(0.81)  # h2o + sub
    plot_args = {"color": "r", "linestyle": "--", "label": "33% coverage"}
    path.plot(**plot_args)
    path.print()

    # for 0.78 coverage
    path = Profile(axes, default_label_color="b")
    path.add_de(0.0)    # co2 + sub
    path.add_de(-0.96)  # co2_fe_bend
    path.add_de(0.31)   # co_form_ts
    path.add_de(-1.03)  # co_fe_o_sub
    path.add_de(1.28)   # co + o_sub
    path.add_de(-0.26)  # h2_fe
    path.add_de(0.61)   # h2o_from1_ts
    path.add_de(-0.37)  # h_fe_h_o
    path.add_de(0.85)   # h2o_form2_ts
    path.add_de(-0.70)  # h2o_sub
    path.add_de(0.97)   # h2o + sub
    plot_args = {"color": "b", "linestyle": "--", "label": "78% coverage"}
    path.plot(**plot_args)
    path.print()

    # Fine adjustments

    # Ticks
    axes.set_xlabel("Reaction coordinate", fontsize="large", weight=font_weight)
    axes.set_ylabel("Energy (kJ/mol)", fontsize="large", weight=font_weight)
    axes.set_xticks([])
    # axes.set_ylim(ymin, ymax)
    axes.tick_params(axis="y", width=axes_width)

    # Spines
    for key in ("top", "bottom", "right"):
        # axes.spines[key].set_visible(False)
        axes.spines[key].set_linewidth(spines_width)
    axes.spines["left"].set_linewidth(spines_width)

    # Legend
    axes.legend(edgecolor="w")

    # Save the figure
    fig.tight_layout()
    fig.savefig("mep.png", dpi=figure_dpi)


if __name__ == "__main__":
    main()
