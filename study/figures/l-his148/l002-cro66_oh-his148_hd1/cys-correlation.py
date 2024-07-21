#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

data_str = "cro65_oh-his146_hd1-dist"
data_label = "Cro66 OH - His148 HD1 Distance [Å]"
corr_str = "ser203_h-asn144_o-dist"
corr_label = r"Cys147 C$_\alpha$ - Cys204 C$_\alpha$ Distance [Å]"


def create_pes(x_data, y_data, plot_title, file_name):
    # Compute 2D histogram
    hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=50)
    hist /= np.sum(hist)  # normalize to probabilities
    hist = np.ma.masked_equal(hist, 0)
    hist = -np.log(hist)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    min_axis = np.floor(min(np.min(x_centers), np.min(y_centers)))
    max_axis = np.ceil(max(np.max(x_centers), np.max(y_centers)))

    # Plot 2D histogram with contours
    plt.contourf(x_centers, y_centers, hist.T, levels=20, cmap="viridis")
    plt.plot(
        [min_axis, max_axis],
        [min_axis, max_axis],
        color="black",
        linestyle="--",
        zorder=1,
    )

    plt.xlim(min_axis, max_axis)
    plt.ylim(min_axis, max_axis)

    plt.colorbar(label="-ln(p)")
    plt.xlabel(corr_label)
    plt.ylabel(data_label)
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Original data
    rogfp_dist_path = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data_str}.npy",
    )
    rogfp_data = np.load(rogfp_dist_path)
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_str}.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_cu_dist_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_str}.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_dist_path)

    # Cys data
    rogfp_cys_data_path = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{corr_str}.npy",
    )
    rogfp_cys_data = np.load(rogfp_cys_data_path)
    rogfp_cys_oxd_data_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{corr_str}.npy",
    )
    rogfp_cys_oxd_data = np.load(rogfp_cys_oxd_data_path)
    rogfp_cys_cu_data_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{corr_str}.npy",
    )
    rogfp_cys_cu_data = np.load(rogfp_cys_cu_data_path)

    create_pes(rogfp_cys_data, rogfp_data, "Reduced", "pes-reduced.png")
    create_pes(rogfp_cys_oxd_data, rogfp_oxd_data, "Oxidized", "pes-oxidized.png")
    create_pes(rogfp_cys_cu_data, rogfp_cu_data, "Cu(I)", "pes-cu.png")
