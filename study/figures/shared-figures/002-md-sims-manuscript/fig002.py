#!/usr/bin/env python3

import os
from xml.etree import ElementTree as ET

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import svgutils.compose as compose
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

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


quantities_to_plot = {
    r"Cys147 C$\alpha$ – Cys204 C$\alpha$": {
        "filename": "cys145_ca-cys202_ca-dist",
        "type": "distance",
        "xlims": (3.0, 7.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
        "ridge": False,
        "label": "A",
        "y_max_factor": 0.1,
    },
    "Gln94 NE2 – Cro66 O3": {
        "filename": "cro65_o3-gln92_ne2-dist",
        "type": "distance",
        "xlims": (2.0, 7.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
        "ridge": False,
        "label": "B",
        "y_max_factor": 0.70,
    },
    r"Cro66 $\Psi$": {
        "filename": "cro65_n3_ca3_c3-val66_n-dihedral",
        "type": "dihedral",
        "xlims": (-120, 240),
        "bin_width": 1.0,
        "bw_method": 0.004,
        "ridge": False,
        "label": "C",
        "y_max_factor": 0.3,
    },
    "Thr203 HG1": {
        "filename": "cro65_oh-thr201_hg1-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
        "ridge": True,
        "label": "D",
    },
    "Tyr145 HH": {
        "filename": "cro65_oh-tyr143_hh-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
        "ridge": True,
        "label": "E",
    },
    "His148 HD1": {
        "filename": "cro65_oh-his146_hd1-dist",
        "type": "distance",
        "xlims": (1.0, 8.0),
        "bin_width": 0.05,
        "bw_method": 0.04,
        "ridge": True,
        "label": "F",
    },
}


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


def add_subfigure_label(ax, label, loc=(0.07, 0.91), fontsize=12, fontweight="normal"):
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


if __name__ == "__main__":
    base_dir = "../../../"

    # Optionally apply a custom style:
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [
        os.path.join(base_dir, "misc/003-figure-style/roboto"),
        os.path.join(base_dir, "misc/003-figure-style/arial"),
    ]
    use_mpl_rc_params(rc_json_path, font_dirs)
    plt.rc("font", family="Arial")
    plt.rc("axes", labelweight="normal")

    # Separate items into ridge vs. other
    ridge_keys = [k for k, v in quantities_to_plot.items() if v.get("ridge", False)]
    other_keys = [k for k, v in quantities_to_plot.items() if not v.get("ridge", False)]

    n_ridge = len(ridge_keys) * 2
    n_other = len(other_keys)
    n_systems = len(paths_system)

    if n_ridge == 0:
        raise ValueError("No 'ridge=True' items found. Please enable at least one.")

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

    # -------------------------------------------------------------------------
    # B. BUILD THE GRIDSPEC
    #
    # We want 2 "major" rows:
    #    - Row 0: a single row for the "other" data (ridge=False).
    #    - Rows 1..n_system: one sub-row per system for the "ridge=True" items.
    #
    # We'll have n_ridge columns in total (for the ridge items).
    # The top row gets the same n_ridge columns, and we distribute "other" items among them.
    # -------------------------------------------------------------------------

    n_buffer = 3
    n_rows = n_systems + n_buffer + 1  # 1 row for "other", plus 1 row per system
    fig = plt.figure(figsize=(7.0, 5.0))
    gs = gridspec.GridSpec(
        nrows=n_rows, ncols=n_ridge, bottom=0.075, top=0.99, left=0.09, right=0.98
    )
    subfigure_label_counter = 0

    # -------------------------------------------------------------------------
    # C. Plot the "other" items (ridge=False) in row=0, distributing them across the n_ridge columns
    # -------------------------------------------------------------------------
    if n_other > 0:
        # Simple approach: each "other" item gets chunk_size columns
        chunk_size = n_ridge // n_other
        leftover = n_ridge % n_other

        col_start = 0
        for i, label_data in enumerate(other_keys):
            col_span = chunk_size + (1 if i < leftover else 0)
            col_end = col_start + col_span
            # Make a subplot that spans [row=0, col_start:col_end]
            ax_other = fig.add_subplot(gs[0, col_start:col_end])

            # Retrieve data
            x_vals = data_all[label_data]["x_vals"]
            pdfs_dict = data_all[label_data]["pdfs"]
            y_max = data_all[label_data]["y_max"]

            # Plot each system in a single subplot
            for sys_lbl in labels_sys_order:
                pdf_y = pdfs_dict[sys_lbl]
                ax_other.plot(
                    x_vals, pdf_y, color=colors_sys[sys_lbl], lw=1.5, label=sys_lbl
                )

            # Aesthetics
            ax_other.set_xlim(*quantities_to_plot[label_data]["xlims"])
            ax_other.set_ylim(
                0, y_max + quantities_to_plot[label_data]["y_max_factor"] * y_max
            )

            # If it's the first "other" item, put a y-axis label
            if i == 0:
                ax_other.set_ylabel("Probability Density")
            ax_other.set_yticks([])

            # Label x-axis depending on distance or dihedral
            q_type = quantities_to_plot[label_data]["type"]
            if q_type == "distance":
                ax_other.set_xlabel(f"{label_data} (Å)")
            else:
                ax_other.set_xlabel(f"{label_data} (degrees)")

            if i == 0:
                ax_other.legend(frameon=False)

            col_start = col_end

            add_subfigure_label(ax_other, quantities_to_plot[label_data]["label"])
            subfigure_label_counter += 1

    # -------------------------------------------------------------------------
    # D. Plot the "ridge=True" items in the rows below (one row per system)
    #
    # Row i_system+1, columns 0..n_ridge-1
    # Each cell is a separate subplot for that (system, ridge_item).
    # We'll do a "filled" style (like your older ridgeline approach).
    # -------------------------------------------------------------------------
    for i_system, sys_lbl in enumerate(labels_sys_order):
        for j, label_data in enumerate(ridge_keys):
            ax_ridge = fig.add_subplot(gs[i_system + n_buffer + 1, j])

            rect = ax_ridge.patch
            rect.set_alpha(0)

            # Retrieve data
            x_vals = data_all[label_data]["x_vals"]
            pdfs_dict = data_all[label_data]["pdfs"]
            y_max = data_all[label_data]["y_max"]

            # This system's PDF
            pdf_y = pdfs_dict[sys_lbl]

            # Plot fill
            ax_ridge.plot(x_vals, pdf_y, color="#373737", lw=0.5)
            ax_ridge.fill_between(
                x_vals, pdf_y, alpha=1, color=colors_sys[sys_lbl], lw=0
            )

            # Aesthetics
            ax_ridge.set_xlim(*quantities_to_plot[label_data]["xlims"])
            ax_ridge.set_ylim(0, y_max)
            ax_ridge.set_yticks([])
            ax_ridge.set_xticks([])

            # Remove spines for a "clean" look
            for spine in ["top", "right", "left", "bottom"]:
                ax_ridge.spines[spine].set_visible(False)

            # If this is the bottom system (the last in labels_sys_order),
            # label the x-axis with the quantity name
            if sys_lbl == labels_sys_order[-1]:
                ax_ridge.set_xticks([2, 4, 6, 8])
                if j == 1:
                    ax_ridge.set_xlabel("Distance from Cro66 (Å)", labelpad=-1)

            # If this is the leftmost column (j==0), we can label the system
            # (like your text() usage or a y-axis label)
            if j == 0:
                ax_ridge.text(
                    quantities_to_plot[label_data]["xlims"][0] - 0.2,
                    0.0,
                    sys_lbl,
                    fontweight="normal",
                    ha="right",
                    va="bottom",
                )

            if i_system == 0:
                add_subfigure_label(
                    ax_ridge,
                    label_data,
                    loc=(0.52, 0.5),
                    fontsize=8,
                    fontweight="normal",
                )
                add_subfigure_label(
                    ax_ridge, quantities_to_plot[label_data]["label"], loc=(0.52, 0.7)
                )
                subfigure_label_counter += 1

    fig.text(0.610, 0.501, "G", fontsize=12, fontweight="normal")

    gs.update(hspace=-0.735)

    plt.savefig("fig002.svg")

    compose.Unit.per_inch["in"] = 1

    compose.Figure(
        "7in",
        "5in",
        compose.SVG("fig002.svg", fix_mpl=True),
        compose.SVG("gfp-relevant-residues.svg")
        .scale(0.65)
        .move(313, 177),  # adjust scale & position
    ).save("fig002.svg")

    tree = ET.parse("fig002.svg")
    root = tree.getroot()

    # Manually set viewBox and fix width/height if needed
    dpi = 600
    root.set("viewBox", f"0 0 {7 * 72} {5 * 72}")
    root.set("width", "7in")
    root.set("height", "5in")

    tree.write("fig002.svg")
