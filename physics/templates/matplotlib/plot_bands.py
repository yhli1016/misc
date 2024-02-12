#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


class Config:
    """Class for holding plotting configurations."""
    def __init__(self) -> None:
        # Figure settings
        self.figure_size = (6.4, 4.8)
        self.figure_dpi = 300
        self.figure_name = "bands.png"

        # Font settings
        self.font_size = 16
        self.font_family = "Liberation Sans"
        self.font_weight = "normal"

        # Axes settings
        self.tick_width = 1.0
        self.tick_length_major = 8
        self.tick_length_minor = 4
        self.spine_width = 1.0

        # Line settings
        self.line_width = 1.5
        self.axline_width = 0.5

    def rc(self, **kwargs) -> None:
        """Change global configurations."""
        plt.rc("font", size=self.font_size, family=self.font_family,
               weight=self.font_weight, **kwargs)

    def figure(self, **kwargs) -> plt.Figure:
        """Create figure from configurations."""
        return plt.figure(figsize=self.figure_size, dpi=self.figure_dpi,
                          **kwargs)


def plot(ax: plt.Axes, config: Config) -> None:
    """Actually plot the data."""
    # Load data
    k_idx = np.load("k_idx.npy")
    k_len = np.load("k_len.npy")
    bands = np.load("bands.npy")

    # Plot data
    num_bands = bands.shape[1]
    for i in range(num_bands):
        ax.plot(k_len, bands[:, i], color="r", linewidth=config.line_width)
    for idx in k_idx:
        ax.axvline(k_len[idx], color="k", linewidth=config.axline_width)

    # Basic ticks settings
    k_label = ["$\Gamma$", "M", "K", "$\Gamma$"]
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


def main():
    config = Config()
    config.rc()

    # Create figure and axes
    fig = config.figure()
    ax = fig.add_subplot()

    # Plot data
    plot(ax, config)

    # Save figure
    fig.tight_layout()
    fig.savefig(config.figure_name)


if __name__ == "__main__":
    main()
