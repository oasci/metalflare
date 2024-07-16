#!/usr/bin/env python3
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

KB = 1.987204259e-3  # kcal/(mol K)
T = 300.0  # Kelvin
bin_min, bin_max = 1, 10
bin_width = 0.2
n_bins = int((bin_max - bin_min) / bin_width)
x_lims = (1.4, 8)
y_lims = (1.4, 8)

fig_label = "b004-pes"
data_x_str = "cro65_oh-his146_hd1-dist"
data_x_label = "Cro66 OH - His148 HD1 Distance [Å]"
data_y_str = "cro65_oh-tyr143_hh-dist"
data_y_label = "Cro66 OH - Tyr145 HH Distance [Å]"


def create_histogram(x_data, y_data, bins=30, in_energy=True):
    hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=bins)
    hist = np.ma.masked_where(hist == 0, hist)
    hist /= np.sum(hist)
    hist = -np.log(hist)
    if in_energy:
        hist *= KB * T
        hist -= np.min(hist)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    bins_info = {"edges": (x_edges, y_edges), "centers": (x_centers, y_centers)}
    return hist, bins_info


def create_pes(x_data, y_data, file_name):
    hist, bins_info = create_histogram(x_data, y_data, bins=n_bins, in_energy=True)

    vmin, vmax = 0, 4.5
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    contour = plt.contourf(
        *bins_info["centers"], hist.T, levels=15, cmap="viridis", norm=norm
    )
    plt.colorbar(
        contour,
        label="PMF [kcal/mol]",
        norm=norm,
        ticks=list(range(vmin, int(np.floor(vmax)) + 1)),
    )

    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


def masked_difference(hist_ref, hist):
    diff = hist - hist_ref

    # Handle the special cases
    mask_not_in_ref = np.where(hist_ref.mask & ~hist.mask)
    mask_only_in_ref = np.where(~hist_ref.mask & hist.mask)
    diff[mask_not_in_ref] = hist[mask_not_in_ref] - np.max(hist)
    diff[mask_only_in_ref] = hist_ref[mask_only_in_ref]

    return diff


def create_difference_pes(data_x_ref, data_y_ref, data_x, data_y, file_name):
    hist_ref, bins_info = create_histogram(
        data_x_ref, data_y_ref, bins=n_bins, in_energy=True
    )
    hist, _ = create_histogram(
        data_x, data_y, bins=bins_info["edges"], in_energy=True
    )

    hist_diff = masked_difference(hist_ref, hist)

    vmin, vmax = -5, 5
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    contour = plt.contourf(
        *bins_info["centers"], hist_diff.T, levels=200, cmap="RdBu_r", norm=norm
    )
    plt.colorbar(
        contour, label="ΔPMF [kcal/mol]", norm=norm, ticks=list(range(vmin, vmax + 1))
    )

    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
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
    path_x_red = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_red = np.load(path_x_red)

    path_x_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_oxd = np.load(path_x_oxd)

    path_x_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_cu = np.load(path_x_cu)

    # Cys data
    path_y_red = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_red = np.load(path_y_red)

    path_y_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_oxd = np.load(path_y_oxd)

    path_y_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_cu = np.load(path_y_cu)

    create_pes(data_x_red, data_y_red, f"{fig_label}-reduced.png")
    create_pes(data_x_oxd, data_y_oxd, f"{fig_label}-oxidized.png")
    create_pes(data_x_cu, data_y_cu, f"{fig_label}-cu.png")

    # Create the difference plot
    create_difference_pes(
        data_x_red, data_y_red, data_x_oxd, data_y_oxd, f"{fig_label}-diff-oxd-red.png"
    )

    create_difference_pes(
        data_x_red, data_y_red, data_x_cu, data_y_cu, f"{fig_label}-diff-cu-red.png"
    )
