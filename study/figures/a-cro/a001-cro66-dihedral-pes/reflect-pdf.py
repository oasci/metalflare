#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params


def load_and_process_data(base_dir, data_str, state):
    path = os.path.join(
        base_dir, f"analysis/{names_state[state]}/data/struct-desc/{data_str}.npy"
    )
    data = np.load(path)
    data = np.degrees(data)  # Convert to degrees
    data = np.concatenate(
        [data, data + 360, data - 360]
    )  # Extend the data for circular KDE
    return data


def circular_kde(data, bw_method=0.001):
    kde = gaussian_kde(data, bw_method=bw_method)
    return kde


def make_dual_circular_pdf_fig(
    x_values,
    kde_dict,
    plt_kwargs,
    x_label="Angle [°]",
    y_label="Density",
    figsize=(3.5, 5.0),
):
    """Create figure with two circular PDFs - one on positive y-axis and one on negative y-axis."""
    fig, ax = plt.subplots(figsize=figsize)

    colors = {
        "Reduced": "#1e2e79",
        "Oxidized": "#EC4067",
        "Cu(I)": "#f99752",
        "Na+": "#1b998b",
    }

    # Plot distributions
    for dtype in kde_dict.keys():
        # Plot x distribution on positive y-axis (wrapped around ±180°)
        kde_values_x = kde_dict[dtype]["x"](x_values)
        ax.plot(
            x_values, kde_values_x, label=f"{dtype}", color=colors[dtype], **plt_kwargs
        )

        # Plot y distribution on negative y-axis (centered around 0°)
        # Create new x values centered around 0
        x_centered = np.linspace(-180, 180, 360)
        kde_values_y = kde_dict[dtype]["y"](x_centered)

        # Reorder the data to center around 0
        mid_point = len(x_centered) // 2
        kde_values_y = np.roll(kde_values_y, mid_point)

        ax.plot(
            x_centered, -kde_values_y, color=colors[dtype], linestyle="--", **plt_kwargs
        )

    # Customize plot
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(-180, 180)

    # Add horizontal line at y=0
    ax.axhline(y=0, color="black", linewidth=0.5)

    # Add vertical gridlines at key angles
    for angle in [-180, -90, 0, 90, 180]:
        ax.axvline(x=angle, color="gray", linestyle=":", alpha=0.3)

    # Adjust legend
    legend = ax.legend(
        frameon=False,
        # bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    # Specify paths and setup
    base_dir = "../../../"

    data_x_str = "cro65_n2_ca2_cb2_cg2-dihedral"
    data_y_str = "cro65_ca2_cb2_cg2_cd1-dihedral"

    names_state = {
        "Reduced": "005-rogfp-glh-md",
        "Oxidized": "007-rogfp-oxd-glh-md",
        "Cu(I)": "006-rogfp-cu-glh-md",
        "Na+": "008-rogfp-na-glh-md",
    }

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Load and process data
    data_types = ["Reduced", "Oxidized", "Cu(I)", "Na+"]
    data_dict = {dtype: {} for dtype in data_types}

    for dtype in data_types:
        data_dict[dtype]["x"] = load_and_process_data(base_dir, data_x_str, dtype)
        data_dict[dtype]["y"] = load_and_process_data(base_dir, data_y_str, dtype)

    # Create KDEs
    kde_dict = {dtype: {} for dtype in data_types}
    for dtype in data_types:
        kde_dict[dtype]["x"] = circular_kde(data_dict[dtype]["x"])
        kde_dict[dtype]["y"] = circular_kde(data_dict[dtype]["y"])

    # Define x values for plotting
    x_values = np.linspace(-180, 180, 360)

    # Create and save the figure
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    fig = make_dual_circular_pdf_fig(
        x_values,
        kde_dict,
        pdf_plt_kwargs,
        x_label="Centered Dihedral Angle [°]",
        y_label="Density",
        figsize=(3.5, 5.0),
    )

    # Add labels for top and bottom distributions
    fig.axes[0].text(
        35, fig.axes[0].get_ylim()[1] * 0.9, "N2-CA2-CB2-CG2", rotation=90, va="top"
    )
    fig.axes[0].text(
        35, fig.axes[0].get_ylim()[0] * 0.9, "CA2-CB2-CG2-CD1", rotation=90, va="bottom"
    )
    fig.axes[0].set_xlim(-40, 40)

    fig.savefig("dual-dihedral-pdf.svg", bbox_inches="tight", dpi=300)
    plt.close()
