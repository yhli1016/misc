#! /usr/bin/env python

from typing import List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from procar import Procar, get_fermi


class Config:
    """Class for holding plotting configurations."""
    def __init__(self) -> None:
        # Figure settings
        self.figure_size = (6.4, 10)
        self.figure_dpi = 300
        self.figure_name = "dos.png"

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


def plot(axes: List[plt.Axes], config: Config) -> None:
    """Actually plot the data."""
    atom_ids = [1]
    symbols = [["x2-y2"], ["dxz", "dyz"]]
    e_step = 0.001
    sigma = 0.05
    strain = ["0.98", "0.99", "1.00", "1.01", "1.02"]
    strain_label = ["-2%", "-1%", "0", "+1%", "+2%"]

    for r in ["rot"]:
        for i, s in enumerate(strain):
            fermi = get_fermi(f"{s}/OUTCAR.{r}")
            e_min = fermi - 1.5
            e_max = fermi + 1.5
            procar = Procar(f"{s}/PROCAR.{r}")
            energy, pdos0 = procar.eval_pdos(atom_ids, symbols[0],
                                             e_min, e_max, e_step, sigma)
            energy, pdos1 = procar.eval_pdos(atom_ids, symbols[1],
                                             e_min, e_max, e_step, sigma)
            energy -= fermi

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


def main():
    config = Config()
    config.rc()

    # Create figure and axes
    fig = config.figure()
    gs = fig.add_gridspec(5, hspace=0)
    axes = gs.subplots(sharex=True, sharey=True)

    # Plot
    plot(axes, config)

    # Save figure
    fig.tight_layout()
    fig.savefig(config.figure_name)


if __name__ == "__main__":
    main()
