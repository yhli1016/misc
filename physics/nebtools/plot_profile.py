#! /usr/bin/env python

from typing import Tuple
from dataclasses import dataclass

import matplotlib.pyplot as plt

from nebtools.profile import EnergyProfile


@dataclass
class Config:
    """Class for holding plotting configurations."""
    # Figure settings
    figure_size: Tuple[float, float] = (7, 5)
    figure_dpi: int = 300
    figure_name: str = "mep_0.png"

    # Font settings
    font_size: int = 12
    font_family: str = "Liberation Sans"
    font_weight: str = "normal"

    # Axes settings
    tick_width: float = 1.5
    tick_length_major: int = 8
    tick_length_minor: int = 4
    spine_width: int = 1.5

    # Line settings
    line_width: float = 0.75
    axline_width: float = 0.5

    def rc(self, **kwargs) -> None:
        """
        Update global configurations.

        There are three approaches to achieve the goal:
            1. call plt.rc
            2. modify plt.rcParams
            3. call plt.rcParams.update
        Anyway, changes to global configurations should be as small as possible.

        :param kwargs: additional global configurations
        :returns: None
        """
        config = {"font.size": self.font_size,
                  "font.family": self.font_family,
                  "font.weight": self.font_weight}
        config.update(kwargs)
        plt.rcParams.update(config)

    def figure(self, **kwargs) -> plt.Figure:
        """
        Create figure from configurations.

        :param kwargs: additional figure configurations
        """
        return plt.figure(figsize=self.figure_size, dpi=self.figure_dpi,
                          **kwargs)


def main():
    # Change global settings and create the figure
    config = Config()
    config.rc()
    fig = config.figure()
    axes = fig.add_subplot()

    # Predefined plotting styles
    # Names for colors: (b)lue, (r)ed, (g)reen, (c)yan, (m)agenta, blac(k), (w)hite
    # Allowed line styles are "-", "--", "-.", ":"
    plot_args = {"linewidth": config.line_width, "linestyle": "--", "color": "k"}
    ref_args = {"linewidth": config.line_width, "linestyle": ":", "color": "k"}
    arrow_args = {"width": 0.001, "head_width": 0.05, "head_length": 5,
                  "length_includes_head": True, "color": "k"}
    text_args = {"color": "k", "size": "x-small"}

    # for X
    path = EnergyProfile(label_color="k", energy_unit="ev")
    path.add_level_delta(0.0)    # co2 + sub
    path.add_level_delta(-1.26)  # co2_fe_bend
    path.add_level_delta(0.39)   # co_form_ts
    path.add_level_delta(-1.23)  # co_fe_o_sub
    path.add_level_delta(1.56)   # co + o_sub
    path.add_level_delta(-0.24)  # h2_fe
    path.add_level_delta(0.41)   # h2o_from1_ts
    path.add_level_delta(-0.31)  # h_fe_h_o
    path.add_level_delta(1.18)   # h2o_form2_ts
    path.add_level_delta(-0.49)  # h2o_sub
    path.add_level_delta(0.69)   # h2o + sub
    plot_args["label"] = "X"
    path.plot(axes, energy_unit="kjm", **plot_args)
    path.add_barrier(axes,2, energy_unit="kjm",
                     ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(axes,3, energy_unit="kjm",
                     ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(axes,6, energy_unit="kjm",
                     ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(axes,8, energy_unit="kjm",
                     ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    print("Fe:")
    path.print(energy_unit="ev")

    # Fine adjustments

    # Ticks
    axes.set_xlabel("Reaction coordinate", fontsize="large",
                    weight=config.font_weight)
    axes.set_ylabel("$\mathrm{Energy (kJ \cdot mol^{-1})}$", fontsize="large",
                    weight=config.font_weight)
    axes.set_xticks([])
    # axes.set_ylim(ymin, ymax)
    axes.tick_params(axis="y", width=config.tick_width)

    # Spines
    for key in ("top", "bottom", "right"):
        # axes.spines[key].set_visible(False)
        axes.spines[key].set_linewidth(config.spine_width)
    axes.spines["left"].set_linewidth(config.spine_width)

    # Legend
    axes.legend(edgecolor="w")

    # Save the figure
    fig.tight_layout()
    fig.savefig(config.figure_name, dpi=config.figure_dpi)


if __name__ == "__main__":
    main()
