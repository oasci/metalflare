#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "g011-pes"
data_x_str = "cro65_og1-glu220_he2-dist"
data_x_label = "Cro66 OG1 - Glu222 HE2 Distance [Å]"
data_y_str = "ser203_og-glu220_he2-dist"
data_y_label = "Glu222 HE2 - Ser205 OG Distance [Å]"

bin_width_angle = 2
bin_width_dist = 0.1
bins = (np.arange(0, 12, bin_width_dist), np.arange(0, 12, bin_width_dist))
x_lims = (1, 9)
x_ticks = np.arange(min(x_lims), max(x_lims) + 0.001, 1)
y_lims = (1, 9)
y_ticks = np.arange(min(y_lims), max(y_lims) + 0.001, 1)


pes_vmin = 0
pes_vmax = 5

pes_diff_vmin = -5
pes_diff_vmax = 5

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

    fig = create_pes(
        data_x_red,
        data_y_red,
        bins=bins,
        vmin=pes_vmin,
        vmax=pes_vmax,
        levels=15,
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
        levels=15,
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
        levels=15,
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
