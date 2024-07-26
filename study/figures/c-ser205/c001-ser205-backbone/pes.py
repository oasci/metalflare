#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes, create_pes_difference

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "c001-pes"
data_x_str = "cys202_c-ser203_n_ca_c-dihedral"
data_x_label = "Ser205 $\phi$ [deg.]"
data_y_str = "ser203_n_ca_c-ala204_n-dihedral"
data_y_label = "Ser205 $\psi$ [deg.]"

bin_width_angle = 2
bin_width_dist = 0.1
bins = (
    np.arange(-360, 0 + 0.01, bin_width_angle),
    np.arange(0, 360 + 0.01, bin_width_angle),
)
x_lims = (-210, -30)
x_ticks = np.arange(-210, 0 + 0.001, 60)
y_lims = (60, 220)
y_ticks = np.arange(min(y_lims), max(y_lims) + 0.001, 60)

pes_vmin = 0
pes_vmax = 3.5

pes_diff_vmin = -4
pes_diff_vmax = 4


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
    data_x_red = np.degrees(data_x_red)
    data_x_red[data_x_red >= 0] -= 360

    path_x_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_oxd = np.load(path_x_oxd)
    data_x_oxd = np.degrees(data_x_oxd)
    data_x_oxd[data_x_oxd >= 0] -= 360

    path_x_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_cu = np.load(path_x_cu)
    data_x_cu = np.degrees(data_x_cu)
    data_x_cu[data_x_cu >= 0] -= 360

    path_y_red = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_red = np.load(path_y_red)
    data_y_red = np.degrees(data_y_red)
    data_y_red[data_y_red <= 0] += 360

    path_y_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_oxd = np.load(path_y_oxd)
    data_y_oxd = np.degrees(data_y_oxd)
    data_y_oxd[data_y_oxd <= 0] += 360

    path_y_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_cu = np.load(path_y_cu)
    data_y_cu = np.degrees(data_y_cu)
    data_y_cu[data_y_cu <= 0] += 360

    fig = create_pes(
        data_x_red,
        data_y_red,
        bins=bins,
        vmin=pes_vmin,
        vmax=pes_vmax,
        levels=13,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-reduced.png")
    plt.close()

    fig = create_pes(
        data_x_oxd,
        data_y_oxd,
        bins=bins,
        vmin=pes_vmin,
        vmax=pes_vmax,
        levels=13,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-oxidized.png")
    plt.close()

    fig = create_pes(
        data_x_cu,
        data_y_cu,
        bins=bins,
        vmin=pes_vmin,
        vmax=pes_vmax,
        levels=13,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-cu.png")
    plt.close()

    fig = create_pes_difference(
        data_x_red,
        data_y_red,
        data_x_oxd,
        data_y_oxd,
        bins=bins,
        vmin=pes_diff_vmin,
        vmax=pes_diff_vmax,
        levels=200,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-diff-oxd-red.png")
    plt.close()

    fig = create_pes_difference(
        data_x_red,
        data_y_red,
        data_x_cu,
        data_y_cu,
        bins=bins,
        vmin=pes_diff_vmin,
        vmax=pes_diff_vmax,
        levels=200,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-diff-cu-red.png")
    plt.close()

    fig = create_pes_difference(
        data_x_oxd,
        data_y_oxd,
        data_x_cu,
        data_y_cu,
        bins=bins,
        vmin=pes_diff_vmin,
        vmax=pes_diff_vmax,
        levels=200,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(*x_lims)
    plt.xticks(x_ticks)
    plt.ylabel(data_y_label)
    plt.ylim(*y_lims)
    plt.yticks(y_ticks)
    plt.tight_layout()
    fig.savefig(f"{fig_label}-diff-cu-oxd.png")
    plt.close()
