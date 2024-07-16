#!/usr/bin/env python3
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "a004-pes"
data1_str = "cro65_cd2_cg2_cb2_ca2-dihedral"
data1_label = "Cro66 CD2-CG2-CB2-CA2 Dihedral [°]"
data2_str = "cro65_cg2_cb2_ca2_c2-dihedral"
data2_label = "Cro66 CG2-CB2-CA2-C2 Dihedral [°]"


def create_histogram(x_data, y_data, bins=50):
    x_bins = np.linspace(-180, 180, bins + 1)
    y_bins = np.linspace(-180, 180, bins + 1)
    hist, _, _ = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins], density=True)
    hist = np.ma.masked_where(hist == 0, hist)
    hist = -np.log(hist)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2
    y_centers = (y_bins[:-1] + y_bins[1:]) / 2
    return hist, x_centers, y_centers


def create_pes(x_data, y_data, file_name):
    hist, x_centers, y_centers = create_histogram(x_data, y_data)

    vmin, vmax = 5, 16
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    contour = ax.contourf(
        x_centers, y_centers, hist.T, levels=100, cmap="viridis", norm=norm
    )
    plt.colorbar(contour, label="-ln(p)", norm=norm, ticks=list(range(vmin, vmax + 1)))

    ax.set_xlabel(data2_label)
    ax.set_xlim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_ylabel(data1_label)
    ax.set_ylim(-180, 180)
    ax.set_yticks(np.arange(-180, 181, 60))
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


def create_difference_pes(x_data_cu, y_data_cu, x_data_red, y_data_red, file_name):
    hist_cu, x_centers, y_centers = create_histogram(x_data_cu, y_data_cu)
    hist_red, _, _ = create_histogram(x_data_red, y_data_red)

    diff_hist = hist_cu - hist_red
    diff_hist = np.ma.masked_where(np.isnan(diff_hist), diff_hist)

    vmin, vmax = -5, 5
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8, 6))
    contour = ax.contourf(
        x_centers, y_centers, diff_hist.T, levels=100, cmap="RdBu_r", norm=norm
    )
    plt.colorbar(contour, label="-Δln(p)", norm=norm, ticks=list(range(vmin, vmax + 1)))

    ax.set_xlabel(data2_label)
    ax.set_xlim(-180, 180)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_ylabel(data1_label)
    ax.set_ylim(-180, 180)
    ax.set_yticks(np.arange(-180, 181, 60))
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


if __name__ == "__main__":
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
    data1_red = np.degrees(data1_red)

    data1_oxd_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data1_str}.npy",
    )
    data1_oxd = np.load(data1_oxd_path)
    data1_oxd = np.degrees(data1_oxd)

    rogfp_cu_dist_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data1_str}.npy",
    )
    data1_cu = np.load(rogfp_cu_dist_path)
    data1_cu = np.degrees(data1_cu)

    # Cys data
    data2_red_path = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_red = np.load(data2_red_path)
    data2_red = np.degrees(data2_red)

    data2_oxd_path = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_oxd = np.load(data2_oxd_path)
    data2_oxd = np.degrees(data2_oxd)

    data2_cu_path = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data2_str}.npy",
    )
    data2_cu = np.load(data2_cu_path)
    data2_cu = np.degrees(data2_cu)

    create_pes(data2_red, data1_red, f"{fig_label}-reduced.png")
    create_pes(data2_oxd, data1_oxd, f"{fig_label}-oxidized.png")
    create_pes(data2_cu, data1_cu, f"{fig_label}-cu.png")

    create_difference_pes(
        data2_cu, data1_cu, data2_red, data1_red, f"{fig_label}-diff-cu-red.png"
    )

    create_difference_pes(
        data2_oxd, data1_oxd, data2_red, data1_red, f"{fig_label}-diff-oxd-red.png"
    )
