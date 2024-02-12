#! /usr/bin/env python

import matplotlib.pyplot as plt

from nebtools.profile import Profile


class Config:
    """Class for holding plotting configurations."""
    def __init__(self) -> None:
        # Figure settings
        self.figure_size = (7, 5)
        self.figure_dpi = 300
        self.figure_name = "mep_0.png"

        # Font settings
        self.font_size = 12
        self.font_family = "Liberation Sans"
        self.font_weight = "normal"

        # Axes settings
        self.tick_width = 1.5
        self.tick_length_major = 8
        self.tick_length_minor = 4
        self.spine_width = 1.5

        # Line settings
        self.line_width = 0.75
        self.axline_width = 0.5

    def rc(self, **kwargs) -> None:
        """Change global configurations."""
        plt.rc("font", size=self.font_size, family=self.font_family,
               weight=self.font_weight, **kwargs)

    def figure(self, **kwargs) -> plt.Figure:
        """Create figure from configurations."""
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
    path = Profile(axes, default_label_color="k", unit="ev")
    path.add_de(0.0)    # co2 + sub
    path.add_de(-1.26)  # co2_fe_bend
    path.add_de(0.39)   # co_form_ts
    path.add_de(-1.23)  # co_fe_o_sub
    path.add_de(1.56)   # co + o_sub
    path.add_de(-0.24)  # h2_fe
    path.add_de(0.41)   # h2o_from1_ts
    path.add_de(-0.31)  # h_fe_h_o
    path.add_de(1.18)   # h2o_form2_ts
    path.add_de(-0.49)  # h2o_sub
    path.add_de(0.69)   # h2o + sub
    plot_args["label"] = "X"
    path.plot(unit="kjm", **plot_args)
    path.add_barrier(2, unit="kjm", ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(3, unit="kjm", ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(6, unit="kjm", ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    path.add_barrier(8, unit="kjm", ref_args=ref_args, arrow_args=arrow_args,
                     text_args=text_args, text_dx=0.1, text_dy=-10)
    print("Fe:")
    path.print(unit="ev")

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
