#! /usr/bin/env python

from typing import List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from plotutils import Config, MultiPlot


@MultiPlot(Config(figure_name="dos.png", figure_size=(6.4, 10)))
def plot(axes: List[plt.Axes], config: Config) -> None:
    """Actually plot the data."""
    strain = ["0.98", "0.99", "1.00", "1.01", "1.02"]
    strain_label = ["-2%", "-1%", "0", "+1%", "+2%"]

    for r in ["rot"]:
        for i, s in enumerate(strain):
            energy = np.load(f"data/pdos/{s}/energy.{r}.npy")
            pdos0 = np.load(f"data/pdos/{s}/pdos0.{r}.npy")
            pdos1 = np.load(f"data/pdos/{s}/pdos1.{r}.npy")

            # Plot
            ax = axes[i]
            ax.plot(energy, pdos0, 'r-', label="$d_{//}$", linewidth=config.line_width)
            ax.plot(energy, pdos1, 'b-', label="$\pi^*$", linewidth=config.line_width)
            ax.axvline(color="k", linestyle="--", linewidth=config.axline_width)
            ax.fill_between(energy, pdos0, where=(energy<=0), color="r", alpha=0.3)
            ax.fill_between(energy, pdos1, where=(energy<=0), color="b", alpha=0.3)
            ax.annotate(f"$\sigma$ = {strain_label[i]}", (0.8, 0.04))

            # Basic ticks settings
            ax.set_xlabel("Energy (eV)", weight=config.font_weight)
            if i == 2:
                ax.set_ylabel("$\mathrm{DOS\ (eV^{-1})}$", weight=config.font_weight)
            ax.set_xlim(np.min(energy), np.max(energy))
            ax.set_ylim(0, 0.05)
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

            # Set legend
            ax.legend(edgecolor="w")


if __name__ == "__main__":
    plot()
