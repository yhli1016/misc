"""Auxiliary functions and classes for matplotlib."""

from typing import Callable, Tuple, List
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


@dataclass
class Config:
    """Class for holding plot configurations."""
    # Figure settings
    figure_size: Tuple[float, float] = (3.2, 2.4)
    figure_dpi: int = 300
    figure_name: str = "fig.png"

    # Font settings
    font_size: int = 8
    font_family: str = "Liberation Sans"
    font_weight: str = "normal"

    # Ticks ans spines settings
    tick_width: float = 0.75
    tick_length_major: float = 4
    tick_length_minor: float = 2
    tick_color: str = "w"
    spine_width: float = 0.75
    spine_color: str = "w"

    # Line settings
    line_width: float = 0.75
    axline_width: float = 0.5
    axline_color: str = "w"

    def rc(self, **kwargs) -> None:
        """
        Update global configurations.

        There are three approaches to achieve the goal:
            1. call plt.rc
            2. modify plt.rcParams
            3. call plt.rcParams.update
        Anyway, changes to global configurations should be as small as possible.

        :param kwargs: additional global configurations
        :returns: None
        """
        config = {"font.size": self.font_size,
                  # "font.family": self.font_family,
                  "font.weight": self.font_weight}
        config.update(kwargs)
        plt.rcParams.update(config)

    def figure(self, **kwargs) -> plt.Figure:
        """
        Create figure from configurations.

        :param kwargs: additional figure configurations
        """
        return plt.figure(figsize=self.figure_size, dpi=self.figure_dpi,
                          **kwargs)


class SinglePlot:
    """
    Decorator for wrapping functions generating a single plot.

    Attributes
    ----------
    config: 'Config' instance
        configurations of the plot
    """
    def __init__(self, config: Config) -> None:
        self.config = config

    def __call__(self, func: Callable) -> Callable:
        def _plot(*args, **kwargs) -> None:
            self.config.rc()
            fig = self.config.figure()
            ax = fig.add_subplot()
            func(ax, self.config, *args, **kwargs)
            fig.tight_layout()
            fig.savefig(self.config.figure_name)
            plt.show()
        return _plot


class MultiPlot(SinglePlot):
    """Decorator for wrapping functions generating a multiple plot."""
    def __call__(self, func: Callable) -> Callable:
        def _plot(*args, **kwargs) -> None:
            self.config.rc()
            fig = self.config.figure()
            gs = fig.add_gridspec(5, hspace=0)
            axes = gs.subplots(sharex=True, sharey=True)
            func(axes, self.config, *args, **kwargs)
            fig.tight_layout()
            fig.savefig(self.config.figure_name)
        return _plot


def plot_single(config: Config, func: Callable, *args, **kwargs) -> None:
    """
    Function form of SinglePlot.

    :param config: configurations of the plot
    :param func: function actually plotting the figure
    :param args: positional arguments for func
    :param kwargs: keyword arguments for func
    :return: None
    """
    config.rc()
    fig = config.figure()
    ax = fig.add_subplot()
    func(ax, config, *args, **kwargs)
    fig.tight_layout()
    fig.savefig(config.figure_name)


def plot_multiple(config: Config, func: Callable, *args, **kwargs) -> None:
    """
    Function form of MultiPlot.

    :param config: configurations of the plot
    :param func: function actually plotting the figure
    :param args: positional arguments for func
    :param kwargs: keyword arguments for func
    :return: None
    """
    config.rc()
    fig = config.figure()
    gs = fig.add_gridspec(5, hspace=0)
    axes = gs.subplots(sharex=True, sharey=True)
    func(axes, config, *args, **kwargs)
    fig.tight_layout()
    fig.savefig(config.figure_name)


@SinglePlot(Config(figure_name="bands_diag.png"))
def plot_diag(ax: plt.Axes,
              config: Config,
              energy: np.ndarray,
              k_len: np.ndarray,
              k_idx: np.ndarray,
              k_label: List[str],
              e_min: float = -1.0,
              e_max: float = 1.0,
              scatter: bool = False) -> None:
    """
    Plot band structure from exact diagonalization.

    :param ax: axes on which the figure will be plot
    :param config: plotting configurations
    :param k_len: (num_kpt,) float64 array
        length of k-path in 1/nm
    :param energy: (num_kpt, num_dos) float64 array
        energies of dos for each k-point in eV
    :param dos: (num_kpt, num_dos) float64 array
        dos in 1/eV
    :param k_len: (num_kpt,) float64 array
        distance of k-path in reciprocal space
    :param k_idx: (num_hsk,) int32 array
        indices of highly-symmetric k-points in k_len
    :param k_label: (num_hsk,) labels of highly-symmetric k-points
    :param e_min: lower bound of energy in plot in eV
    :param e_max: upper bound of energy in plot in eV
    :param scatter: whether to plot band data as scatter plot
    :return: None
    """
    # Filter data
    idx = np.where((energy[0] >= e_min) & (energy[0] <= e_max))[0]
    idx_min, idx_max = idx[0], idx[-1]
    energy = energy[:, idx_min:idx_max]

    # Plot data
    num_bands = energy.shape[1]
    for i in range(num_bands):
        if scatter:
            ax.scatter(k_len, energy[:, i], color="r", s=0.2)
        else:
            ax.plot(k_len, energy[:, i], color="r",
                    linewidth=config.line_width)
    for idx in k_idx:
        ax.axvline(k_len[idx], color="k", linewidth=config.axline_width)

    # Basic ticks settings
    # ax.set_xlabel()
    ax.set_ylabel("Energy (eV)")
    ax.set_xlim(0, np.amax(k_len))
    ax.set_ylim(e_min, e_max)
    ax.set_xticks(k_len[k_idx])
    # ax.set_yticks()
    ax.set_xticklabels(k_label)
    # ax.set_yticklabels()

    # Advanced ticks settings
    ax.minorticks_on()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator())
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator())
    # ax.yaxis.set_major_locator(ticker.MultipleLocator())
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(which="both", direction="in", width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)
    ax.tick_params(axis="x", which="both", length=0)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)


@SinglePlot(Config(figure_name="bands_tbpm.png"))
def plot_tbpm(ax: plt.Axes,
              config: Config,
              energy: np.ndarray,
              dos: np.ndarray,
              k_len: np.ndarray,
              k_idx: np.ndarray,
              k_label: List[str],
              e_min: float = -1.0,
              e_max: float = 1.0,
              num_grid: Tuple[int, int] = (200, 200),
              cmap: str = "viridis",
              **kwargs) -> None:
    """
    Plot band structure from TBPM DOS.

    :param ax: axes on which the figure will be plot
    :param config: plotting configurations
    :param k_len: (num_kpt,) float64 array
        length of k-path in 1/nm
    :param energy: (num_kpt, num_dos) float64 array
        energies of dos for each k-point in eV
    :param dos: (num_kpt, num_dos) float64 array
        dos in 1/eV
    :param k_len: (num_kpt,) float64 array
        distance of k-path in reciprocal space
    :param k_idx: (num_hsk,) int32 array
        indices of highly-symmetric k-points in k_len
    :param k_label: (num_hsk,) labels of highly-symmetric k-points
    :param e_min: lower bound of energy in plot in eV
    :param e_max: upper bound of energy in plot in eV
    :param num_grid: (nx, ny) number of grid-points for interpolation along
        x and y directions when plotting the wave function
    :param cmap: color map for plotting the wave function
    :param kwargs: parameters for plt.scatter() and plt.imshow()
    :return: None
    """
    # Convert arrays to vectors for interpolation
    idx = np.where((energy[0] >= e_min) & (energy[0] <= e_max))[0]
    idx_min, idx_max = idx[0], idx[-1]
    energy = energy[:, idx_min:idx_max]
    dos = dos[:, idx_min:idx_max]
    x = k_len.repeat(energy.shape[1])
    y = energy.flatten()
    z = dos.flatten()

    # Interpolate and plot data
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    x_fi = np.linspace(x_min, x_max, num_grid[0])
    y_fi = np.linspace(y_min, y_max, num_grid[1])
    x_grid, y_grid = np.meshgrid(x_fi, y_fi)
    xy_fi = np.c_[x_grid.ravel(), y_grid.ravel()]
    z_fi = griddata((x, y), z, xy_fi, method="cubic")
    z_fi = z_fi.reshape(num_grid)
    extent = (x_min, x_max, y_min, y_max)
    ax.imshow(z_fi, cmap=cmap, interpolation="none", origin="lower",
              extent=extent, **kwargs)
    for idx in k_idx:
        ax.axvline(k_len[idx], color=config.axline_color,
                   linewidth=config.axline_width)

    # Adjustment
    ax.set_xticks(k_len[k_idx])
    ax.set_xticklabels(k_label)
    ax.set_ylabel("Energy (eV)")
    ax.set_aspect((x_max - x_min)/(y_max - y_min))
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", color=config.tick_color,
                   width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)
    ax.tick_params(axis="x", which="both", length=0)
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)
        ax.spines[pos].set_color(config.spine_color)
