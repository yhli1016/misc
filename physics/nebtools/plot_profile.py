#! /usr/bin/env python

import matplotlib.pyplot as plt

from nebtools.profile import Profile


def main():
    # Figure Settings
    figure_size = (7, 5)
    figure_dpi = 300

    # Font settings
    font_size = 12
    font_family = "Arial"
    # font_family = "Liberation Sans"
    font_weight = "normal"

    # Axes and spines settings
    # ymin = -2.0
    # ymax = 2.0
    axes_width = 1.5
    spines_width = 1.5

    # Line settings
    line_width = 0.75

    # Change global settings and create the figure
    plt.rc("font", size=font_size, family=font_family, weight=font_weight)
    fig, axes = plt.subplots(figsize=figure_size)

    # Names for colors: (b)lue, (r)ed, (g)reen, (c)yan, (m)agenta, blac(k), (w)hite
    # Allowed line styles are "-", "--", "-.", ":"

    # for X
    path = Profile(axes, default_label_color="b", unit="ev")
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
    plot_args = {"linewidth": line_width, "linestyle": "--", "color": "b",
                 "label": "X"}
    path.plot(unit="kjm", **plot_args)
    print("Fe:")
    path.print(unit="ev")

    # Fine adjustments

    # Ticks
    axes.set_xlabel("Reaction coordinate", fontsize="large", weight=font_weight)
    axes.set_ylabel("$\mathrm{Energy (kJ \cdot mol^{-1})}$", fontsize="large", weight=font_weight)
    axes.set_xticks([])
    # axes.set_ylim(ymin, ymax)
    axes.tick_params(axis="y", width=axes_width)

    # Spines
    for key in ("top", "bottom", "right"):
        # axes.spines[key].set_visible(False)
        axes.spines[key].set_linewidth(spines_width)
    axes.spines["left"].set_linewidth(spines_width)

    # Legend
    axes.legend(edgecolor="w")

    # Save the figure
    fig.tight_layout()
    fig.savefig("mep_0.png", dpi=figure_dpi)


if __name__ == "__main__":
    main()
