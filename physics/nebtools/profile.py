"""Class for evaluating and plotting energy profiles."""

from abc import ABC, abstractmethod

import numpy as np
from matplotlib.collections import LineCollection


class Path:
    """
    Class for evaluating the energy profile of a single reaction path.

    Attributes
    ----------
    label: List[string]
        labels for reactants and products
    energy: List[float]
        energies of reactants and products
    unit: string
        unit of energy
    """
    def __init__(self, unit="ev") -> None:
        """
        :param unit: string
            unit of energy, should be either "ev" or "kjm"
        :raises ValueError: if unit is neither "ev" nor "kjm"
        """
        self.label = []
        self.energy = []
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        self.unit = unit

    def scale_energy(self, unit="ev"):
        """
        Get scaled energy in given unit.

        :param unit: string
            unit of energy, should be either "ev" or "kjm"
        """
        ev2kjm = 96.4916
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        if unit == "kjm":
            if self.unit == "kjm":
                scale_factor = 1.0
            else:
                scale_factor = ev2kjm
        else:
            if self.unit == "ev":
                scale_factor = 1.0
            else:
                scale_factor = 1.0 / ev2kjm
        energy = np.array(self.energy) * scale_factor
        return energy

    def add_eng(self, label, energy):
        """
        Add an energy level in the reaction path.

        :param string label: label for the state
        :param float energy: energy of the state
        :return: None
        """
        if energy is not None:  # DO NOT DELETE THIS LINE!
            self.label.append(label)
            self.energy.append(energy)

    def eval_eng(self, unit="eV"):
        """
        Print energy levels and differences of the reaction path.

        :param string unit: unit of energies for output
        :return: None
        """
        energy = self.scale_energy(unit=unit)
        for i, label in enumerate(self.label):
            eng = energy[i]
            eng_align = eng - energy[0]
            eng_delta = eng - energy[i-1] if i > 0 else 0
            if i > 0:
                print("%16s : %8.2f%8.2f" % (label, eng_align, eng_delta))
            else:
                print("%16s : %8.2f%8s" % (label, eng_align, "diff"))


class MultiPath(ABC):
    """
    Base class for evaluating the energy profile of multiple reaction paths.

    Attributes
    ----------
    paths: List[Path]
        list of single reaction paths
    """
    def __init__(self) -> None:
        self.paths = []

    @abstractmethod
    def gen_paths(self):
        """To be implemented in derived classes."""
        pass

    def eval_eng(self, unit="ev"):
        """
        Print energy levels and differences of the reaction paths.

        :param string unit: unit of energies for output
        :return: None
        """
        self.gen_paths()
        for path in self.paths:
            path.eval_eng(unit=unit)
            print()


class Profile:
    """
    Class for plotting the energy profile of a reaction path.

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
    unit: string
        unit of energy
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
                 label_length=0.6, label_width=1.8, unit="ev") -> None:
        self.axes = axes
        self.default_label_color = default_label_color
        self.default_connector = default_connector
        self.label_length = label_length
        self.label_width = label_width
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        self.unit = unit
        self.react_coord = []
        self.energy = []
        self.label_color = []
        self.connector = []

    @property
    def curr_x(self):
        """Get the reaction coordinate of last energy level."""
        if len(self.react_coord) == 0:
            return -1
        else:
            return self.react_coord[-1]

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

    def scale_energy(self, unit="ev"):
        """
        Get scaled energy in given unit.

        :param unit: string
            unit of energy, should be either "ev" or "kjm"
        """
        ev2kjm = 96.4916
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        if unit == "kjm":
            if self.unit == "kjm":
                scale_factor = 1.0
            else:
                scale_factor = ev2kjm
        else:
            if self.unit == "ev":
                scale_factor = 1.0
            else:
                scale_factor = 1.0 / ev2kjm
        energy = np.array(self.energy) * scale_factor
        return energy

    def plot(self, unit="ev", **kwargs):
        """
        Plot energy profile onto axes.

        :param str unit: unit for energy levels
        :param dict kwargs: arguments to be passed to Axes.plot()
        :return: None
        """
        # Prepare data
        react_coord = np.array(self.react_coord)
        energy = self.scale_energy(unit=unit)

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
        self.axes.plot(react_coord_fi, energy_fi, **kwargs)

        # Add energy levels
        levels = []
        for i, coord in enumerate(react_coord):
            p0 = (coord - 0.5 * self.label_length, energy[i])
            p1 = (coord + 0.5 * self.label_length, energy[i])
            levels.append((p0, p1))
        levels = LineCollection(levels, colors=self.label_color,
                                linewidth=self.label_width)
        self.axes.add_collection(levels)

    def add_barrier(self, i=1, unit="ev", ref_args=None,
                    arrow_dx=0.0, arrow_args=None,
                    text_dx=0.0, text_dy=0.0, text_args=None):
        """
        Add a barrier notation to the energy profile.

        :param int i: index of energy level for which the notation will
            be added with respect to the previous energy level
        :param str unit: unit for energy levels
        :param dict ref_args: settings of the reference line
        :param float arrow_dx: shift of arrow from x1
        :param dict arrow_args: settings of the arrow
        :param float text_dx: shift of text from the center of arrow along x
        :param float text_dy: shift of text from the center of arrow along y
        :param dict text_args: settings of the text
        :return: None
        """
        # Prepare data
        react_coord = np.array(self.react_coord)
        energy = self.scale_energy(unit=unit)

        # Draw the horizontal reference line
        x0, x1 = react_coord[i-1], react_coord[i]
        y0, y1 = energy[i-1], energy[i]
        p0 = x0 + 0.5 * self.label_length
        p1 = x1 + 0.5 * self.label_length
        self.axes.plot((p0, p1), (y0, y0), **ref_args)

        # Draw the vertical arrow
        ax, ay = x1 + arrow_dx, y0
        dx, dy = 0, y1 - y0
        self.axes.arrow(ax, ay, dx, dy, **arrow_args)

        # Add the text
        tx = ax + text_dx
        ty = y0 + dy * 0.5 + text_dy
        self.axes.text(tx, ty, f"{dy:.2f}", **text_args)

    def print(self, unit="ev"):
        """
        Print absolute energy levels and energy differences of the profile.

        :param str unit: unit for the energy levels
        :return: None
        """
        energy = self.scale_energy(unit=unit)
        for i, coord in enumerate(self.react_coord):
            eng_delta = energy[i] - energy[i-1] if i > 0 else 0
            if i > 0:
                print("%4d : %8.2f%8.2f" % (coord, energy[i], eng_delta))
            else:
                print("%4d : %8.2f%8s" % (coord, energy[i], "diff"))
