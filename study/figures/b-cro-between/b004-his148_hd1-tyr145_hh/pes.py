#!/usr/bin/env python3
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "b004-pes"
data1_str = "cro65_oh-his146_hd1-dist"
data1_label = "Cro66 OH - His148 HD1 Distance [Å]"
data2_str = "cro65_oh-tyr143_hh-dist"
data2_label = "Cro66 OH - Tyr145 HH Distance [Å]"


def create_histogram(x_data, y_data, bins=50):
    hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=bins, density=True)
    hist = np.ma.masked_where(hist == 0, hist)
    hist = -np.log(hist)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    return hist, x_centers, y_centers


def create_pes(x_data, y_data, plot_title, file_name):
    hist, x_centers, y_centers = create_histogram(x_data, y_data)

    vmin, vmax = 0, 8
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    contour = plt.contourf(
        x_centers, y_centers, hist.T, levels=100, cmap="viridis", norm=norm
    )
    plt.colorbar(contour, label="-ln(p)", norm=norm, ticks=list(range(vmin, vmax + 1)))

    plt.xlabel(data2_label)
    plt.xlim(1.5, 6.5)
    plt.ylabel(data1_label)
    plt.ylim(1.5, 6.5)
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


def create_difference_pes(x_data_cu, y_data_cu, x_data_red, y_data_red, file_name):
    hist_cu, x_centers, y_centers = create_histogram(x_data_cu, y_data_cu)
    hist_red, _, _ = create_histogram(x_data_red, y_data_red)

    # Fill masked values with a high number (less probable than any real data point)
    high_value = np.max(hist_cu.filled(0)) + 1
    hist_cu_filled = hist_cu.filled(high_value)
    hist_red_filled = hist_red.filled(high_value)

    diff_hist = hist_cu_filled - hist_red_filled

    # Mask the difference where both original histograms were masked
    diff_hist = np.ma.masked_where((hist_cu.mask) & (hist_red.mask), diff_hist)

    vmin, vmax = -5, 5
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    contour = plt.contourf(
        x_centers, y_centers, diff_hist.T, levels=100, cmap="RdBu_r", norm=norm
    )
    plt.colorbar(contour, label="-Δln(p)", norm=norm, ticks=list(range(vmin, vmax + 1)))

    plt.xlabel(data2_label)
    plt.xlim(1.5, 6.5)
    plt.ylabel(data1_label)
    plt.ylim(1.5, 6.5)
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
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data1_str}.npy",
    )
    data1_red = np.load(rogfp_dist_path)

    data1_oxd_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data1_str}.npy",
    )
    data1_oxd = np.load(data1_oxd_path)

    rogfp_cu_dist_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data1_str}.npy",
    )
    data1_cu = np.load(rogfp_cu_dist_path)

    # Cys data
    data2_red_path = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_red = np.load(data2_red_path)

    data2_oxd_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_oxd = np.load(data2_oxd_path)

    data2_cu_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_cu = np.load(data2_cu_path)

    create_pes(data2_red, data1_red, "Reduced", f"{fig_label}-reduced.png")
    create_pes(data2_oxd, data1_oxd, "Oxidized", f"{fig_label}-oxidized.png")
    create_pes(data2_cu, data1_cu, "Cu(I)", f"{fig_label}-cu.png")

    # Create the difference plot (Cu - Reduced)
    create_difference_pes(
        data2_cu, data1_cu, data2_red, data1_red, f"{fig_label}-diff-cu-red.png"
    )

    create_difference_pes(
        data2_oxd, data1_oxd, data2_red, data1_red, f"{fig_label}-diff-oxd-red.png"
    )
