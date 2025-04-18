#!/usr/bin/env python3

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# If you want to use your custom style:
# from metalflare.analysis.figures import use_mpl_rc_params

###############################################################################
# 1) DEFINE YOUR PLOT CONFIGURATIONS
###############################################################################

# System paths and their color mapping
paths_system = {
    "Reduced": "005-rogfp-glh-md",
    "Oxidized": "007-rogfp-oxd-glh-md",
    "Cu(I)": "006-rogfp-cu-glh-md",
}
labels_sys_order = ["Reduced", "Oxidized", "Cu(I)"]
colors_sys = {
    "Reduced": "#1e2e79",
    "Oxidized": "#EC4067",
    "Cu(I)": "#f99752",
}

# ----------------------------------------------------------------
# The main dictionary describing *all* quantities to plot:
#  - "filename" is the basename of your .npy file (no .npy extension).
#  - "type" can be "distance" or "dihedral".
#  - "xlims" is your desired plot range.
#  - Optional "bin_width" and "bw_method" for the KDE can be set individually.
# ----------------------------------------------------------------

quantities_to_plot = {
    "Gln94 NE2 - Cro66 O3": {
        "filename": "cro65_o3-gln92_ne2-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
    },
    "Thr64 O - Cro66 N1": {
        "filename": "cro65_n1-thr62_o-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
    },
    "Glu222 HE2 - Cro66 OG1": {
        "filename": "cro65_og1-glu220_he2-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
    },
}

###############################################################################
# 2) HELPER FUNCTIONS
###############################################################################


def load_raw_data(base_dir, path_sys, q_info):
    """
    Loads raw data from a .npy file. If 'type' is 'dihedral', it handles
    conversion from radians to degrees (if needed) and duplication by ±360.
    """
    file_path = os.path.join(
        base_dir, f"analysis/{path_sys}/data/struct-desc/{q_info['filename']}.npy"
    )
    data = np.load(file_path)

    if q_info["type"] == "dihedral":
        # If your angles are already in degrees, remove np.degrees
        data = np.degrees(data)
        # Duplicate ±360 for better continuity
        data = np.concatenate([data, data + 360, data - 360])

    return data


def compute_kde_for_quantity(all_data, q_info):
    """
    Given all raw data from multiple systems for one quantity,
    compute x-values (linspace) and the PDFs via Gaussian KDE.
    """
    x_min, x_max = q_info["xlims"]
    bin_width = q_info.get("bin_width", 0.05)
    bw_method = q_info.get("bw_method", 0.04)

    # Create a linear space for the range
    n_bins = int((x_max - x_min) / bin_width)
    x_values = np.linspace(x_min, x_max, n_bins)

    results = {}
    y_max = 0.0

    for sys_label, data_array in all_data.items():
        kde = gaussian_kde(data_array, bw_method=bw_method)
        y_values = kde(x_values)
        results[sys_label] = y_values
        y_here = np.max(y_values)
        if y_here > y_max:
            y_max = y_here

    return x_values, results, y_max


def add_subfigure_label(ax, label, loc=(0.075, 0.95), fontsize=12, fontweight="bold"):
    """Add subfigure label (e.g., A, B, C) to a given axis."""
    ax.text(
        *loc,
        label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight=fontweight,
        va="center",
        ha="center",
    )


###############################################################################
# 3) MAIN SCRIPT
###############################################################################

if __name__ == "__main__":
    base_dir = "../../../"

    # Optionally apply a custom style:
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto"), os.path.join(base_dir, "misc/003-figure-style/arial")]
    use_mpl_rc_params(rc_json_path, font_dirs)
    plt.rc("font", family="Arial")

    # Separate items into ridge vs. other
    label_keys = [k for k, v in quantities_to_plot.items()]
    n_systems = len(paths_system)

    # -------------------------------------------------------------------------
    # A. Load all data, compute KDE for each quantity
    # -------------------------------------------------------------------------
    data_all = {}
    for q_label, q_info in quantities_to_plot.items():
        # Collect raw data for each system
        per_system_data = {}
        for sys_lbl, sys_path in paths_system.items():
            arr = load_raw_data(base_dir, sys_path, q_info)
            per_system_data[sys_lbl] = arr

        x_vals, pdfs, y_max = compute_kde_for_quantity(per_system_data, q_info)
        data_all[q_label] = {
            "x_vals": x_vals,
            "pdfs": pdfs,
            "y_max": y_max,
        }


    fig, axes = plt.subplots(nrows=1, ncols=len(label_keys), figsize=(4.5, 3.0))
    subfigure_label_counter = 0



    col_start = 0
    for i, label_data in enumerate(label_keys):
        # Make a subplot that spans [row=0, col_start:col_end]
        ax = axes[i]

        # Retrieve data
        x_vals = data_all[label_data]["x_vals"]
        pdfs_dict = data_all[label_data]["pdfs"]
        y_max = data_all[label_data]["y_max"]

        # Plot each system in a single subplot
        for sys_lbl in labels_sys_order:
            pdf_y = pdfs_dict[sys_lbl]
            ax.plot(
                x_vals, pdf_y, color=colors_sys[sys_lbl], lw=1.5, label=sys_lbl
            )

        # Aesthetics
        ax.set_xlim(*quantities_to_plot[label_data]["xlims"])
        ax.set_xticks([2.0, 4.0, 6.0])
        ax.set_ylim(0, y_max + 0.1 * y_max)

        # If it's the first "other" item, put a y-axis label
        if i == 0:
            ax.set_ylabel("Probability Density")
        ax.set_yticks([])

        # Label x-axis depending on distance or dihedral
        ax.set_xlabel(f"{label_data} [Å]", fontsize=8)

        if i == len(label_keys)-1:
            ax.legend(frameon=False)

        col_start += 1

        add_subfigure_label(ax, chr(65 + subfigure_label_counter))
        subfigure_label_counter += 1

    fig.tight_layout()
    plt.savefig("fig007.svg")

