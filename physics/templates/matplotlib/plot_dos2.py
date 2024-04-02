#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from plotutils import Config, plot_single


def plot(ax: plt.Axes, config: Config) -> None:
    """
    Plot DOS as function of energy.

    :param ax: axes on which DOS will be plotted
    :param config: configurations of the plot
    :returns: None
    """
    data_dir = "data/dos"
    energy = np.load(f"{data_dir}/energy_dos.npy")
    dos = np.load(f"{data_dir}/dos.npy")
    ax.plot(energy, dos, 'r-', linewidth=config.line_width)
    # ax.scatter(energy, dos, 0.5, 'r')

    # Basic ticks settings
    ax.set_xlabel("Energy (eV)", fontsize="large", weight=config.font_weight)
    ax.set_ylabel("DOS ($eV^{-1}$)", fontsize="large", weight=config.font_weight)
    # ax.set_xlim(-12, 7)
    # ax.set_ylim(0, 0.18)
    # ax.set_xicks()
    # ax.set_yticks()
    # ax.set_xticklabels()
    # ax.set_yticklabels()

    # Advanced ticks settings
    ax.minorticks_on()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
    ax.tick_params(which="both", direction="in", width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)
    # ax.grid(linestyle='--', linewidth=config.line_width*0.5)


if __name__ == "__main__":
    plot_single(Config(figure_name="dos2.png"), plot)
