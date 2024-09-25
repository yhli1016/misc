"""Class for evaluating and plotting energy profiles."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


EV2KJM = 96.4916


@dataclass
class ReactionState:
    """
    Class representing a single initial, transition or final state in the
    reaction.
    """
    label: str = "X"
    energy: float = 0.0


class ReactionPath:
    """
    Class representing a single reaction path consisting of multiple states.

    Attributes
    ----------
    states: List[ReactionState]
        list of states on the reaction path
    energy_unit: str
        unit of energy
    """
    def __init__(self, energy_unit: str = "ev") -> None:
        """
        :param energy_unit: unit of energy, should be either "ev" or "kjm"
        """
        assert energy_unit in ("kjm", "ev")
        self.states = []
        self.energy_unit = energy_unit

    def get_scaled_energies(self, energy_unit: str = "ev") -> np.ndarray:
        """
        Get scaled energies in given unit.

        :param energy_unit: unit of energy, should be either "ev" or "kjm"
        :return: scaled energies as an array
        """
        assert energy_unit in ("kjm", "ev")
        if energy_unit == "kjm" and self.energy_unit == "ev":
            scale_factor = EV2KJM
        elif energy_unit == "ev" and self.energy_unit == "kjm":
            scale_factor = 1.0 / EV2KJM
        else:
            scale_factor = 1.0
        energies = np.array([_.energy for _ in self.states]) * scale_factor
        return energies

    def add_state(self, label: str, energy: float) -> None:
        """
        Add an energy level in the reaction path.

        :param label: label for the state
        :param energy: energy of the state
        :return: None
        """
        # When used in PathCollection, energies passed to this function may be
        # None due to missing or uninitialized data. In that case, drop this
        # state.
        if energy is not None:  # DO NOT DELETE THIS LINE!
            self.states.append(ReactionState(label, energy))

    def eval_energy_differences(self, energy_unit: str = "eV") -> None:
        """
        Print energy levels and differences of the reaction path.

        :param energy_unit: unit of energies for output
        :return: None
        """
        energy = self.get_scaled_energies(energy_unit=energy_unit)
        for i, state in enumerate(self.states):
            label = state.label
            eng_align = energy[i] - energy[0]
            eng_delta = energy[i] - energy[i-1] if i > 0 else 0
            if i > 0:
                print(f"{label:>16s} : {eng_align:8.2f}{eng_delta:8.2f}")
            else:
                print(f"{label:>16s} : {eng_align:8.2f}{'diff':>8s}")


class PathCollection(ABC):
    """
    Base class for evaluating the energy differences of multiple reaction paths.

    The user must implement the build_paths method in their dericed class.

    Attributes
    ----------
    paths: List[ReactionPath]
        list of reaction paths
    """
    def __init__(self) -> None:
        self.paths = []

    def eval_energy_differences(self, energy_unit: str = "ev") -> None:
        """
        Print energy levels and differences of the reaction paths.

        :param energy_unit: unit of energies for output
        :return: None
        """
        self.build_paths()
        for path in self.paths:
            path.eval_energy_differences(energy_unit=energy_unit)
            print()

    @abstractmethod
    def build_paths(self) -> None:
        """Should be implemented in derived classes."""
        pass


@dataclass
class EnergyLevel:
    """
    Class representing a single energy level in the energy profile.
    """
    position: float = 0.0
    energy: float = 0.0
    label_length: float = 0.6  # not used for now
    label_width: float = 1.8  # not used for now
    label_color: str = "k" # not used for now
    connector: int = 1


@dataclass
class EnergyProfile:
    """
    Class for plotting the energy profile of a reaction path.

    Attributes
    ---------
    energy_levels: List[EnergyLevel]
        energy levels in the file
    energy_unit: str
        unit of energies in energy_levels
    label_length: float
        default length of labels of energy levels
    label_width: float
        default width of labels of energy levels
    label_color: str
        default color of labels
    connector: integer
        default order of polynomial for connecting labels (energy levels)
        in the profile
    """
    def __init__(self,
                 energy_unit: str = "ev",
                 label_length: float = 0.6,
                 label_width: float = 1.8,
                 label_color: str = "k",
                 connector: int = 1) -> None:
        assert energy_unit in ("kjm", "ev")
        self.energy_levels = []
        self.energy_unit = energy_unit
        self.label_length = label_length
        self.label_width = label_width
        self.label_color = label_color
        self.connector = connector

    @property
    def curr_pos(self) -> float:
        """Get the reaction coordinate of last energy level."""
        if len(self.energy_levels) == 0:
            return -1.0
        else:
            return self.energy_levels[-1].position

    def add_level(self,
                  position: float,
                  energy: float,
                  label_length: float = None,
                  label_width: float = None,
                  label_color: str = None,
                  connector: int = None) -> None:
        """
        Add an absolute energy level to the profile.

        :param position: x-coordinate of the energy level
        :param energy: y-coordinate of the energy level
        :param label_length: length of label
        :param label_width: width of label
        :param label_color: color of the energy level
        :param connector: order of polynomial for the label
        :return: None
        """
        if label_length is None:
            label_length = self.label_length
        if label_width is None:
            label_width = self.label_width
        if label_color is None:
            label_color = self.label_color
        if connector is None:
            connector = self.connector
        level = EnergyLevel(position, energy,
                            label_length, label_width, label_color,
                            connector)
        self.energy_levels.append(level)

    def add_level_delta(self, dy: float, dx: float = 1.0, **kwargs: dict) -> None:
        """
        Add a relative energy level to the profile.

        :param dy: relative energy w.r.t. previous energy level
        :param dx: incremental of reaction coordinate of the energy level
        :param kwargs: arguments to be passed to the 'add_eng' method
        :return: None
        """
        position = self.curr_pos + dx
        if len(self.energy_levels) == 0:
            energy = dy
        else:
            energy = self.energy_levels[-1].energy + dy
        self.add_level(position, energy, **kwargs)

    def linear_fit(self,
                   xco: np.ndarray,
                   yco: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Generate straight connectors."""
        x0 = xco[0] + 0.5 * self.label_length
        x1 = xco[1] - 0.5 * self.label_length
        y0, y1 = yco[0], yco[1]
        xfi = np.linspace(x0, x1, 51)
        yfi = np.linspace(y0, y1, 51)
        return xfi, yfi

    @staticmethod
    def cubic_fit(xco: np.ndarray,
                  yco: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Generate cubic connectors."""
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

    def get_scaled_energies(self, energy_unit: str = "ev") -> np.ndarray:
        """
        Get scaled energies in given unit.

        :param energy_unit: unit of energy, should be either "ev" or "kjm"
        :return: scaled energies as an array
        """
        assert energy_unit in ("kjm", "ev")
        if energy_unit == "kjm" and self.energy_unit == "ev":
            scale_factor = EV2KJM
        elif energy_unit == "ev" and self.energy_unit == "kjm":
            scale_factor = 1.0 / EV2KJM
        else:
            scale_factor = 1.0
        energies = np.array([_.energy for _ in self.energy_levels]) * scale_factor
        return energies

    def plot(self, axes: plt.axes, energy_unit: str = "ev", **kwargs: dict) -> None:
        """
        Plot energy profile onto axes.

        :param axes: axes onto which the profile will be plotted
        :param energy_unit: unit for energy levels
        :param kwargs: arguments to be passed to Axes.plot()
        :return: None
        """
        # Prepare data
        position = np.array([_.position for _ in self.energy_levels])
        energy = self.get_scaled_energies(energy_unit)

        # Draw the connectors
        pos_fi = np.array([])
        energy_fi = np.array([])
        for i in range(energy.size-1):
            if self.energy_levels[i].connector == 1:
                fit_func = self.linear_fit
            else:
                fit_func = self.cubic_fit
            xfi, yfi = fit_func(position[i:i+2], energy[i:i+2])
            pos_fi = np.append(pos_fi, xfi)
            energy_fi = np.append(energy_fi, yfi)
        axes.plot(pos_fi, energy_fi, **kwargs)

        # Add energy levels
        levels = []
        for i in range(position.shape[0]):
            pos = position.item(i)
            p0 = (pos - 0.5 * self.label_length, energy.item(i))
            p1 = (pos + 0.5 * self.label_length, energy.item(i))
            levels.append((p0, p1))
        levels = LineCollection(levels, colors=self.label_color,
                                linewidth=self.label_width)
        axes.add_collection(levels)

    def add_barrier(self,
                    axes: plt.Axes,
                    i: int = 1,
                    energy_unit: str = "ev",
                    ref_args: dict = None,
                    arrow_dx: float = 0.0,
                    arrow_args: dict = None,
                    text_dx: float = 0.0,
                    text_dy: float = 0.0,
                    text_args: dict = None) -> None:
        """
        Add a barrier notation to the energy profile.

        :param axes: axes onto which the profile will be plotted
        :param i: index of energy level for which the notation will
            be added with respect to the previous energy level
        :param energy_unit: unit for energy levels
        :param ref_args: settings of the reference line
        :param float arrow_dx: shift of arrow from x1
        :param arrow_args: settings of the arrow
        :param text_dx: shift of text from the center of arrow along x
        :param text_dy: shift of text from the center of arrow along y
        :param text_args: settings of the text
        :return: None
        """
        # Prepare data
        position = np.array([_.position for _ in self.energy_levels])
        energy = self.get_scaled_energies(energy_unit)

        # Draw the horizontal reference line
        x0, x1 = position.item(i-1), position.item(i)
        y0, y1 = energy.item(i-1), energy.item(i)
        p0 = x0 + 0.5 * self.label_length
        p1 = x1 + 0.5 * self.label_length
        axes.plot((p0, p1), (y0, y0), **ref_args)

        # Draw the vertical arrow
        ax, ay = x1 + arrow_dx, y0
        dx, dy = 0, y1 - y0
        axes.arrow(ax, ay, dx, dy, **arrow_args)

        # Add the text
        tx = ax + text_dx
        ty = y0 + dy * 0.5 + text_dy
        axes.text(tx, ty, f"{dy:.2f}", **text_args)

    def print(self, energy_unit: str = "ev") -> None:
        """
        Print absolute energy levels and energy differences of the profile.

        :param energy_unit: unit for the energy levels
        :return: None
        """
        position = np.array([_.position for _ in self.energy_levels])
        energy = self.get_scaled_energies(energy_unit)
        for i in range(position.shape[0]):
            pos = int(position[i])
            eng_delta = energy[i] - energy[i-1] if i > 0 else 0
            if i > 0:
                print(f"{pos:4d} : {energy[i]:8.2f}{eng_delta:8.2f}")
            else:
                print(f"{pos:4d} : {energy[i]:8.2f}{'diff':>8s}")
