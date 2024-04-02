#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from plotutils import Config, SinglePlot


@SinglePlot(Config(figure_name="bands.png"))
def plot(ax: plt.Axes, config: Config) -> None:
    """Actually plot the data."""
    # Load data
    data_dir = "data/bands"
    k_idx = np.load(f"{data_dir}/k_idx.npy")
    k_len = np.load(f"{data_dir}/k_len.npy")
    bands = np.load(f"{data_dir}/bands.npy")

    # Plot data
    num_bands = bands.shape[1]
    for i in range(num_bands):
        ax.plot(k_len, bands[:, i], color="r", linewidth=config.line_width)
    for idx in k_idx:
        ax.axvline(k_len[idx], color="k", linewidth=config.axline_width)

    # Basic ticks settings
    k_label = ["L", "G", "X", "W", "K", "G"]
    # ax.set_xlabel()
    ax.set_ylabel("Energy (eV)")
    ax.set_xlim(0, np.amax(k_len))
    # ax.set_ylim()
    ax.set_xticks(k_len[k_idx])
    # ax.set_yticks()
    ax.set_xticklabels(k_label)
    # ax.set_yticklabels()

    # Advanced ticks settings
    # ax.minorticks_on()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator())
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator())
    # ax.yaxis.set_major_locator(ticker.MultipleLocator())
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(which="both", direction="in", width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)
    ax.tick_params(axis="x", which="both", length=0)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)


if __name__ == "__main__":
    plot()
