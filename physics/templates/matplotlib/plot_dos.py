#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def load_data():
    x = np.load("energy.npy")
    y = np.load("dos.npy")
    return x, y


def main():
    # Figure settings
    figure_size = (6.4, 4.8)
    figure_dpi = 300
    figure_name = "dos.png"

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
    energy, dos = load_data()
    ax.plot(energy, dos, 'r-', linewidth=line_width)

    # Basic ticks settings
    ax.set_xlabel("Energy (eV)", fontsize="large", weight=font_weight)
    ax.set_ylabel("DOS (1/eV)", fontsize="large", weight=font_weight)
    ax.set_xlim(-9, 9)
    ax.set_ylim(0, 0.18)
    # ax.set_xicks()
    # ax.set_yticks()
    # ax.set_xticklabels()
    # ax.set_yticklabels()

    # Advanced ticks settings
    ax.minorticks_on()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
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
