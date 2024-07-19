#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import numpy as np
from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes, create_pes_difference

os.chdir(os.path.dirname(os.path.realpath(__file__)))

bin_min, bin_max = 1, 10
bin_width = 0.2
n_bins = int((bin_max - bin_min) / bin_width)
x_lims = (3, 7)
x_ticks = np.arange(2, 7 + 1, 1)
y_lims = (1.4, 8)
y_ticks = np.arange(2, 8 + 1, 1)

fig_label = "f008-pes"
data_x_str = "cys145_ca-cys202_ca-dist"
data_x_label = r"Cys147 C$_\alpha$ - Cys204 C$_\alpha$ Distance [Å]"
data_y_str = "cys145_sg-cys202_sg-dist"
data_y_label = r"Cys147 SG - Cys204 SG Distance [Å]"


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
        data_x_red, data_y_red, bins=n_bins, vmin=0, vmax=4.5, levels=15, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(3, 6)
    plt.xticks(np.arange(3, 6 + 1, 1))
    plt.ylabel(data_y_label)
    plt.ylim(3.2, 8.2)
    plt.yticks(np.arange(3, 8 + 1, 1))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-reduced.png")
    plt.close()

    fig = create_pes(
        data_x_oxd, data_y_oxd, bins=n_bins, vmin=0, vmax=4.5, levels=15, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(3, 7)
    plt.xticks(np.arange(3, 7 + 0.001, 1))
    plt.ylabel(data_y_label)
    plt.ylim(1, 5)
    plt.yticks(np.arange(1, 5 + 0.001, 1))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-oxidized.png")
    plt.close()

    fig = create_pes(
        data_x_cu, data_y_cu, bins=n_bins, vmin=0, vmax=4.5, levels=15, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(3, 7)
    plt.xticks(np.arange(3, 7 + 0.001, 1))
    plt.ylabel(data_y_label)
    plt.ylim(1, 5)
    plt.yticks(np.arange(1, 5 + 0.001, 1))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-cu.png")
    plt.close()

    bins = np.arange(1, 8, bin_width)
    fig = create_pes_difference(
        data_x_oxd,
        data_y_oxd,
        data_x_cu,
        data_y_cu,
        bins=bins,
        vmin=-6,
        vmax=6,
        levels=200,
        T=300.0,
    )
    plt.xlabel(data_x_label)
    plt.xlim(3, 7)
    plt.xticks(np.arange(3, 7 + 0.001, 1))
    plt.ylabel(data_y_label)
    plt.ylim(1, 5)
    plt.yticks(np.arange(1, 5 + 0.001, 1))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-diff-cu-oxd.png")
    plt.close()
