#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def load_data():
    k_idx = np.load("k_idx.npy")
    k_len = np.load("k_len.npy")
    bands = np.load("bands.npy")
    return k_idx, k_len, bands


def main():
    # Figure settings
    figure_size = (6.4, 4.8)
    figure_dpi = 300
    figure_name = "bands.png"

    # Font settings
    font_size = 16
    font_family = "Arial"
    # font_family = "Liberation Sans"
    font_weight = "normal"

    # Line settings
    line_width = 1.5
    tick_width = 1.0
    tick_length_major = 8
    tick_length_minor = 4

    # Change global settings
    plt.rc("font", size=font_size, family=font_family, weight=font_weight)

    # Create figure and axes
    fig = plt.figure(figsize=figure_size, dpi=figure_dpi)
    ax = fig.add_subplot()

    # Plot data
    k_idx, k_len, bands = load_data()
    num_bands = bands.shape[1]
    for i in range(num_bands):
        ax.plot(k_len, bands[:, i], color="r", linewidth=line_width)
    for idx in k_idx:
        ax.axvline(k_len[idx], color="k", linewidth=tick_width)

    # Basic ticks settings
    k_label = ["$\Gamma$", "M", "K", "$\Gamma$"]
    # ax.set_xlabel()
    ax.set_ylabel("Energy (eV)")
    ax.set_xlim((0, np.amax(k_len)))
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
    ax.tick_params(which="both", direction="in", width=tick_width)
    ax.tick_params(which="major", length=tick_length_major)
    ax.tick_params(which="minor", length=tick_length_minor)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(tick_width)

    # Save figure
    fig.tight_layout()
    fig.savefig(figure_name)


if __name__ == "__main__":
    main()
