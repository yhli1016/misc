#! /usr/bin/env python
"""Plot carrier mobility without decorator or wrapper."""

import numpy as np
import matplotlib.pyplot as plt

from plotutils import Config


def index_eng(energies, eng):
    return np.where(energies <= eng)[0][-1]


def plot_dos_mob(energy_dos: np.ndarray, dos: np.ndarray,
                 energy_mob: np.ndarray, mobility: np.ndarray,
                 mobility_cutoff: float = 400) -> None:
    """
    Plot DOS and carrier density as functions of energy in a single graph.

    :param energy_dos: energies of DOS in eV
    :param dos: DOS in 1/eV
    :param energy_mob: energies of mobility in eV
    :param mobility: carrier mobility in cm^2/(Vs)
    :param mobility_cutoff: lower bound for plotting mobility in cm^2/(Vs)
    :returns: None
    """
    # Create config, figure and axes
    config = Config(figure_name="dos_mob3.png")
    fig = config.figure()
    ax1 = fig.add_subplot()

    # Mobility w.r.t. energy
    mob_idx = np.where(mobility>=mobility_cutoff)[0]
    eng_plot, mob_plot = energy_mob[mob_idx], mobility[mob_idx]
    ax1.scatter(eng_plot, mob_plot, s=3, c="r")
    ax1.set_xlabel("Energy (eV)", fontsize="large", weight=config.font_weight)
    ax1.set_ylabel("Mobility ($V^{-1}s^{-1}cm^2$)", fontsize="large",
                   weight=config.font_weight, color="r")
    ax1.set_ylim(np.min(mob_plot)*0.95, np.max(mob_plot)*1.05)

    # DOs w.r.t. energy
    ax2 = ax1.twinx()
    dos_idx = np.where((energy_dos>=energy_mob[0])&(energy_dos<=energy_mob[-1]))
    eng_plot, dos_plot = energy_dos[dos_idx], dos[dos_idx]
    ax2.plot(energy_dos[dos_idx], dos[dos_idx], 'b--', linewidth=config.line_width)
    ax2.set_ylabel("DOS ($eV^{-1}$)", fontsize="large", weight=config.font_weight,
                   color="b")
    ax2.set_ylim(0, np.max(dos_plot)*1.05)

    # Advanced ticks settings
    ax1.spines["left"].set_color("r")
    ax1.tick_params(axis="y", which="both", colors="r")
    ax2.spines["left"].set_color("r")  # This has to be done for ax2, otherwise it is still black.
    ax2.spines["right"].set_color("b")
    ax2.tick_params(axis="y", which="both", colors="b")
    for ax in (ax1, ax2):
        ax.minorticks_on()
        ax.tick_params(which="both", direction="in", width=config.tick_width)
        ax.tick_params(which="major", length=config.tick_length_major)
        ax.tick_params(which="minor", length=config.tick_length_minor)

        # Set spines
        for pos in ("top", "bottom", "left", "right"):
            ax.spines[pos].set_linewidth(config.spine_width)

    # Show and save
    plt.show()
    fig.tight_layout()
    fig.savefig(config.figure_name)


def main():
    data_dir = "data/mobility"
    energy_dos = np.load(f"{data_dir}/energy_dos.npy")
    dos = np.load(f"{data_dir}/dos.npy")
    energy_dc = np.load(f"{data_dir}/energy_dc.npy")
    mobility = np.load(f"{data_dir}/mobility.npy")
    plot_dos_mob(energy_dos, dos, energy_dc, mobility)


if __name__ == "__main__":
    main()
