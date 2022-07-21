#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def load_data():
    pass


def main():
    # Figure settings
    figure_size = (6.4, 4.8)
    figure_dpi = 300
    figure_name = "fig.png"

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
    x, y = load_data()
    ax.plot(x, y, 'r-', linewidth=line_width)

    # Basic ticks settings
    # ax.set_xlabel("", fontsize="large", weight=font_weight)
    # ax.set_ylabel("", fontsize="large", weight=font_weight)
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
