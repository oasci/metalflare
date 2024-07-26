#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes, create_pes_difference

os.chdir(os.path.dirname(os.path.realpath(__file__)))

bin_min, bin_max = -180, 180
bin_width = 3
n_bins = int((bin_max - bin_min) / bin_width)
bins = (np.linspace(-180, 180, n_bins + 1), np.linspace(-180, 180, n_bins + 1))

fig_label = "a001-pes"
data_x_str = "cro65_n2_ca2_cb2_cg2-dihedral"
data_x_label = "Cro66 N2-CA2-CB2-CG2 Dihedral [°]"
data_y_str = "cro65_ca2_cb2_cg2_cd1-dihedral"
data_y_label = "Cro66 CA2-CB2-CG2-CD1 Dihedral [°]"


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

    path_x_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_oxd = np.load(path_x_oxd)
    data_x_oxd = np.degrees(data_x_oxd)

    path_x_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_x_str}.npy",
    )
    data_x_cu = np.load(path_x_cu)
    data_x_cu = np.degrees(data_x_cu)

    path_y_red = os.path.join(
        base_dir,
        f"analysis/005-rogfp-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_red = np.load(path_y_red)
    data_y_red = np.degrees(data_y_red)
    data_y_red[data_y_red < 0] += 360
    data_y_red[data_y_red > 0] -= 180

    path_y_oxd = os.path.join(
        base_dir,
        f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_oxd = np.load(path_y_oxd)
    data_y_oxd = np.degrees(data_y_oxd)
    data_y_oxd[data_y_oxd < 0] += 360
    data_y_oxd[data_y_oxd > 0] -= 180

    path_y_cu = os.path.join(
        base_dir,
        f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_y_str}.npy",
    )
    data_y_cu = np.load(path_y_cu)
    data_y_cu = np.degrees(data_y_cu)
    data_y_cu[data_y_cu < 0] += 360
    data_y_cu[data_y_cu > 0] -= 180

    fig = create_pes(
        data_x_red, data_y_red, bins=bins, vmin=0, vmax=5, levels=13, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(-60, 60)
    plt.xticks(np.arange(-60, 61, 30))
    plt.ylabel(data_y_label)
    plt.ylim(-60, 60)
    plt.yticks(np.arange(-60, 61, 30))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-reduced.png")
    plt.close()

    fig = create_pes(
        data_x_oxd, data_y_oxd, bins=bins, vmin=0, vmax=5, levels=13, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(-60, 60)
    plt.xticks(np.arange(-60, 61, 30))
    plt.ylabel(data_y_label)
    plt.ylim(-60, 60)
    plt.yticks(np.arange(-60, 61, 30))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-oxidized.png")
    plt.close()

    fig = create_pes(
        data_x_cu, data_y_cu, bins=bins, vmin=0, vmax=5, levels=13, T=300.0
    )
    plt.xlabel(data_x_label)
    plt.xlim(-60, 60)
    plt.xticks(np.arange(-60, 61, 30))
    plt.ylabel(data_y_label)
    plt.ylim(-60, 60)
    plt.yticks(np.arange(-60, 61, 30))
    plt.tight_layout()
    fig.savefig(f"{fig_label}-cu.png")
    plt.close()

    # fig = create_pes_difference(
    #     data_x_red,
    #     data_y_red,
    #     data_x_oxd,
    #     data_y_oxd,
    #     bins=bins,
    #     vmin=-5,
    #     vmax=5,
    #     levels=200,
    #     T=300.0,
    # )
    # plt.xlabel(data_x_label)
    # plt.xlim(-60, 60)
    # plt.xticks(np.arange(-60, 61, 30))
    # plt.ylabel(data_y_label)
    # plt.ylim(-60, 60)
    # plt.yticks(np.arange(-60, 61, 30))
    # plt.tight_layout()
    # fig.savefig(f"{fig_label}-diff-oxd-red.png")
    # plt.close()

    # fig = create_pes_difference(
    #     data_x_red,
    #     data_y_red,
    #     data_x_cu,
    #     data_y_cu,
    #     bins=bins,
    #     vmin=-5,
    #     vmax=5,
    #     levels=200,
    #     T=300.0,
    # )
    # plt.xlabel(data_x_label)
    # plt.xlim(-60, 60)
    # plt.xticks(np.arange(-60, 61, 30))
    # plt.ylabel(data_y_label)
    # plt.ylim(-60, 60)
    # plt.yticks(np.arange(-60, 61, 30))
    # plt.tight_layout()
    # fig.savefig(f"{fig_label}-diff-cu-red.png")
    # plt.close()
