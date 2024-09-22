#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes

os.chdir(os.path.dirname(os.path.realpath(__file__)))

bin_min, bin_max = 1, 10
bin_width = 0.1
n_bins = int((bin_max - bin_min) / bin_width)
x_lims = (1.4, 6)
x_ticks = np.arange(2, 6 + 0.01, 1)
y_lims = (1.4, 7)
y_ticks = np.arange(2, 7 + 0.01, 1)

fig_label = "g001-pes-combined"
data_x_str = "cro65_oh-his146_hd1-dist"
data_x_label = "His148 HD1 to Cro66 OH Distance [Å]"
data_y_str = "cro65_oh-tyr143_hh-dist"
data_y_label = "Tyr145 HH to Cro66 OH Distance [Å]"

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Load the data
    paths_x = [
        os.path.join(
            base_dir, f"analysis/005-rogfp-glh-md/data/struct-desc/{data_x_str}.npy"
        ),
        os.path.join(
            base_dir, f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_x_str}.npy"
        ),
        os.path.join(
            base_dir, f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_x_str}.npy"
        ),
    ]
    paths_y = [
        os.path.join(
            base_dir, f"analysis/005-rogfp-glh-md/data/struct-desc/{data_y_str}.npy"
        ),
        os.path.join(
            base_dir, f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_y_str}.npy"
        ),
        os.path.join(
            base_dir, f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_y_str}.npy"
        ),
    ]

    data_x = [np.load(path) for path in paths_x]
    data_y = [np.load(path) for path in paths_y]

    fig, axes = plt.subplots(
        3, 1, figsize=(3.5, 8), sharex=True, sharey=True
    )  # Share the x-axis
    labels = ["A", "B", "C"]

    titles = ["Reduced", "Oxidized", "Cu"]
    for i, (ax, x_data, y_data, label) in enumerate(zip(axes, data_x, data_y, labels)):
        fig = create_pes(
            x_data,
            y_data,
            bins=n_bins,
            vmin=0,
            vmax=4,
            levels=15,
            T=300.0,
            ax=ax,
            colorbar=False,
        )
        ax.set_xlim(*x_lims)
        ax.set_xticks(x_ticks)
        ax.set_ylim(*y_lims)
        ax.set_yticks(y_ticks)
        ax.tick_params(axis="both", which="major", labelsize=8)
        ax.text(
            0.98,
            0.97,
            label,
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
        )

        # if i == 1:  # Only display the colorbar for the middle plot
        #     cbar = fig.colorbar(ax.collections[0], ax=ax, orientation='vertical')
        #     cbar.set_label("PMF [kcal/mol]")

    # Set the shared labels
    fig.text(0.55, 0.01, data_x_label, ha="center", fontsize=8, fontweight="bold")
    fig.text(
        0.02,
        0.52,
        data_y_label,
        va="center",
        rotation="vertical",
        fontsize=8,
        fontweight="bold",
    )

    plt.tight_layout(rect=[0.05, 0.01, 1.0, 1.0])
    fig.savefig(f"{fig_label}.png", dpi=1000)  # Save the combined figure
    plt.close()
