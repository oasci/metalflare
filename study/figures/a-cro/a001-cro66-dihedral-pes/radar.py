#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

bin_min, bin_max = -180, 180
bin_width = 3
n_bins = int((bin_max - bin_min) / bin_width)
bins = (np.linspace(-180, 180, n_bins + 1), np.linspace(-180, 180, n_bins + 1))

data_x_str = "cro65_n2_ca2_cb2_cg2-dihedral"
data_x_label = "Cro66 N2-CA2-CB2-CG2 Dihedral [°]"
data_y_str = "cro65_ca2_cb2_cg2_cd1-dihedral"
data_y_label = "Cro66 CA2-CB2-CG2-CD1 Dihedral [°]"

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
    "na": "008-rogfp-na-glh-md",
}


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


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Load and process data
    data_types = ["reduced", "oxidized", "cu", "na"]
    data_dict = {dtype: {} for dtype in data_types}

    for dtype in data_types:
        data_dict[dtype]["x"] = load_and_process_data(base_dir, data_x_str, dtype)
        data_dict[dtype]["y"] = load_and_process_data(base_dir, data_y_str, dtype)

    # Create KDEs
    kde_dict = {dtype: {} for dtype in data_types}
    for dtype in data_types:
        kde_dict[dtype]["x"] = circular_kde(data_dict[dtype]["x"])
        kde_dict[dtype]["y"] = circular_kde(data_dict[dtype]["y"])

    # Plotting
    fig, ax = plt.subplots(figsize=(3.5, 3.5), subplot_kw=dict(projection="polar"))
    colors = {
        "reduced": "#1e2e79",
        "oxidized": "#EC4067",
        "cu": "#f99752",
        "na": "#1b998b",
    }
    angles = np.linspace(0, 2 * np.pi, 360)  # Go from 0 to 2π

    for dtype in data_types:
        kde_values_x = kde_dict[dtype]["x"](np.degrees(angles))
        kde_values_y = kde_dict[dtype]["y"](np.degrees(angles))
        ax.plot(
            angles,
            kde_values_x,
            label=f"{dtype.capitalize()} - {data_x_label}",
            color=colors[dtype],
            linewidth=1.5,
        )
        ax.plot(
            angles,
            kde_values_y,
            label=f"{dtype.capitalize()} - {data_y_label}",
            color=colors[dtype],
            linewidth=1.5,
            linestyle="--",
        )

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, ax.get_ylim()[1])  # Set lower limit to 0
    ax.set_yticks([])  # Remove radial ticks
    ax.set_xticks(np.arange(0, 2 * np.pi, 1 / 2 * np.pi))  # Go from 0 to 2π
    ax.set_xticklabels(["0°", "90°", "180°", "270°"])

    plt.tight_layout()
    plt.savefig("cro66-dihedral-circular.svg", bbox_inches="tight", dpi=300)
    plt.close()
