#! /usr/bin/env python

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from plotutils import Config, plot_single


def plot(ax: plt.Axes, config: Config) -> None:
    """Actually plot the data."""
    data_dir = "data/ebs"
    # Load KLABELS
    hsp = np.loadtxt(f"{data_dir}/KLABELS", dtype=np.string_, skiprows=1,
                     usecols=(0, 1))
    labels_x = [float(i) for i in hsp[:-1, 1].tolist()]
    group_labels = hsp[:-1, 0].tolist()
    group_labels = [i.decode('utf-8', 'ignore') for i in group_labels]
    for i, label in enumerate(group_labels):
        if re.search(r"GAMMA", label) is not None:
            group_labels[i] = r"$\Gamma$"
        elif re.search(r"\w+_\d+", label) is not None:
            group_labels[i] = r"$\mathrm{%s}$" % label

    # Load band data
    data = np.loadtxt(f"{data_dir}/EBS.dat")
    kpt = data[:, 0]
    eng = data[:, 1]
    wgt = data[:, 2]

    # Plot
    img = ax.scatter(kpt, eng, c=wgt, s=wgt*2, cmap="YlOrRd")
    for x in labels_x[1:-1]:
        ax.axvline(x, color="k", linewidth=config.axline_width)
    ax.axhline(0.0, color="k", linewidth=config.axline_width, linestyle="--")
    plt.gcf().colorbar(img)

    # Basic ticks settings
    # ax.set_xlabel()
    ax.set_ylabel("Energy (eV)")
    ax.set_xlim(0, np.amax(kpt))
    ax.set_ylim(-2, 2)
    ax.set_xticks(labels_x)
    # ax.set_yticks()
    ax.set_xticklabels(group_labels)
    # ax.set_yticklabels()

    # Advanced ticks settings
    # ax.minorticks_on()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator())
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator())
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.tick_params(which="both", direction="in", width=config.tick_width)
    ax.tick_params(which="major", length=config.tick_length_major)
    ax.tick_params(which="minor", length=config.tick_length_minor)
    ax.tick_params(axis="x", which="both", length=0)

    # Set spines
    for pos in ("top", "bottom", "left", "right"):
        ax.spines[pos].set_linewidth(config.spine_width)


if __name__ == "__main__":
    plot_single(Config(figure_name="ebs2.png"), plot)
