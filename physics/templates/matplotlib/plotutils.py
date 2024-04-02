"""Auxiliary functions and classes for matplotlib."""

from typing import Callable, Tuple
from dataclasses import dataclass

import matplotlib.pyplot as plt


@dataclass
class Config:
    """Class for holding plot configurations."""
    # Figure settings
    figure_size: Tuple[float, float] = (6.4, 4.8)
    figure_dpi: int = 300
    figure_name: str = "fig.png"

    # Font settings
    font_size: int = 16
    font_family: str = "Liberation Sans"
    font_weight: str = "normal"

    # Ticks ans spines settings
    tick_width: float = 1.0
    tick_length_major: float = 8
    tick_length_minor: float = 4
    spine_width: float = 1.0

    # Line settings
    line_width: float = 1.5
    axline_width: float = 0.5

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
                  "font.family": self.font_family,
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
            func(ax, self.config, *args, *kwargs)
            fig.tight_layout()
            fig.savefig(self.config.figure_name)
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
    func(ax, config, *args, *kwargs)
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
