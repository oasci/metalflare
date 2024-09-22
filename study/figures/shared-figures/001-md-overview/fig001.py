#!/usr/bin/env python3

import os

import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pdfs import compute_pdf

os.chdir(os.path.dirname(os.path.realpath(__file__)))


def load_data(base_dir, file_paths):
    """Load numpy data from the given file paths."""
    return [np.load(os.path.join(base_dir, path)) for path in file_paths]


def compute_pdfs(data_list, x_values, bw_method):
    """Compute PDFs for each dataset in data_list."""
    return [compute_pdf(data, x_values, bw_method=bw_method) for data in data_list]


def plot_pdf(
    ax,
    x_values,
    pdfs,
    labels,
    colors,
    x_label,
    y_label,
    plot_x_bounds,
    plot_y_bounds,
    linewidth=1.5,
):
    """Plot PDFs on a given axis."""
    for pdf, label, color in zip(pdfs, labels, colors):
        ax.plot(
            x_values,
            pdf,
            label=label,
            color=color,
            alpha=1.0,
            linewidth=linewidth,
            linestyle="-",
        )
    ax.set_xlim(plot_x_bounds)
    ax.set_ylim(plot_y_bounds)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def insert_img(ax, img_path):
    """Insert an SVG image into the specified axis."""
    img = image.imread(img_path)
    ax.imshow(img)
    ax.axis("off")


def create_main_figure():
    """Create the main figure with four subplots."""
    fig, axs = plt.subplots(2, 2, figsize=(6, 4.5))
    return fig, axs


def add_legend(ax, colors, labels):
    """Add a legend to the given axis."""
    handles = [plt.Line2D([0], [0], color=color, lw=2) for color in colors]
    ax.legend(handles, labels, frameon=False)


def main():
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Define data file paths for each PDF
    cys_paths = [
        "analysis/005-rogfp-glh-md/data/struct-desc/cys145_ca-cys202_ca-dist.npy",
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cys145_ca-cys202_ca-dist.npy",
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cys145_ca-cys202_ca-dist.npy",
        "analysis/008-rogfp-na-glh-md/data/struct-desc/cys145_ca-cys202_ca-dist.npy",
    ]
    thr_paths = [
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
        "analysis/008-rogfp-na-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
    ]
    ser_glu_paths = [
        "analysis/005-rogfp-glh-md/data/struct-desc/ser203_og-glu220_he2-dist.npy",
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/ser203_og-glu220_he2-dist.npy",
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/ser203_og-glu220_he2-dist.npy",
        "analysis/008-rogfp-na-glh-md/data/struct-desc/ser203_og-glu220_he2-dist.npy",
    ]

    # Load data
    cys_data = load_data(base_dir, cys_paths)
    thr_data = load_data(base_dir, thr_paths)
    ser_glu_data = load_data(base_dir, ser_glu_paths)

    # Compute PDFs
    x_bounds = (1, 10)
    bin_width = 0.05
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.05  # Manually tuned

    cys_pdfs = compute_pdfs(cys_data, x_values, bw_method)
    thr_pdfs = compute_pdfs(thr_data, x_values, bw_method)
    ser_glu_pdfs = compute_pdfs(ser_glu_data, x_values, bw_method)

    # Colors and labels
    colors = ["#1E2E79", "#EC4067", "#f99752", "#D3C4E3"]  # Example color scheme
    labels = ["Reduced", "Oxidized", "Cu(I)", "Na+"]

    # Create the main figure and subplots
    fig, axs = create_main_figure()

    # Plot Cys PDF (top-right)
    plot_pdf(
        axs[0, 1],
        x_values,
        cys_pdfs,
        labels,
        colors,
        r"Cys147 C$_\alpha$ - Cys204 C$_\alpha$ Distance [Å]",
        "Probability density",
        (3, 7),
        (0, None),
    )

    # Plot Thr PDF (bottom-left)
    plot_pdf(
        axs[1, 0],
        x_values,
        thr_pdfs,
        labels,
        colors,
        "Cro66 OH - Thr203 HG1 Distance [Å]",
        "Probability density",
        (1, 7.5),
        (0, None),
    )

    # Plot Ser-Glu PDF (bottom-right)
    plot_pdf(
        axs[1, 1],
        x_values,
        ser_glu_pdfs,
        labels,
        colors,
        "Glu222 HE2 - Ser205 OG Distance [Å]",
        "Probability density",
        (1, 8),
        (0, None),
    )

    # Add legend only in the top-left panel
    add_legend(axs[0, 1], colors, labels)

    fig.tight_layout()

    # Add sim setup
    sim_img_path = "sim.png"
    insert_img(axs[0, 0], sim_img_path)
    axs[0, 0].set_position([0.005, 0.52, 0.55, 0.46])  # [left, bottom, width, height]

    # Add subplot labels
    left_column_x = 0.12
    right_column_x = 0.61

    top_row_y = 0.95
    buttom_row_y = 0.46

    fig.text(
        left_column_x,
        top_row_y,
        "A",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color="black",
    )
    fig.text(
        right_column_x,
        top_row_y,
        "B",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color="black",
    )
    fig.text(
        left_column_x,
        buttom_row_y,
        "C",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color="black",
    )
    fig.text(
        right_column_x,
        buttom_row_y,
        "D",
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
        color="black",
    )

    # Save the figure
    fig.savefig("fig001.svg")
    fig.savefig("fig001.png", dpi=300)


if __name__ == "__main__":
    main()
