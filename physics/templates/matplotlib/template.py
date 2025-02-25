#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from plotutils import Config, SinglePlot


@SinglePlot(Config(figure_name="fig.png"))
def plot(ax: plt.Axes, config: Config) -> None:
    """Actually plot the data."""
    x, y = np.load("x.npy"), np.load("y.npy")
    ax.plot(x, y, 'r-', linewidth=config.line_width)

    # Basic ticks settings
    # ax.set_xlabel("", fontsize="large", weight=config.font_weight)
    # ax.set_ylabel("", fontsize="large", weight=config.font_weight)
    # ax.set_xlim()
    # ax.set_ylim()
    # ax.set_xicks()
    # ax.set_yticks()
    # ax.set_xticklabels()
    # ax.set_yticklabels()

    # Advanced ticks settings
    # ax.minorticks_on()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator())
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator())
    # ax.yaxis.set_major_locator(ticker.MultipleLocator())
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator())
    ax.tick_params(which="both", direction="in", width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)


if __name__ == "__main__":
    plot()
