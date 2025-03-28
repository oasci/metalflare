#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as grid_spec
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))


paths_data = {
    # "Glu222 HE2": "cro65_og1-glu220_he2-dist",
    # "Thr203 HG1": "cro65_oh-thr201_hg1-dist",
    "His148 HD1": "cro65_oh-his146_hd1-dist",
    "Tyr145 HH": "cro65_oh-tyr143_hh-dist",
    "Thr64 O": "cro65_n1-thr62_o-dist",
    "Gln94 NE2": "cro65_o3-gln92_ne2-dist",
}
x_lims_data = {
    "Glu222 HE2": (1.0, 8.0),
    "Thr203 HG1": (1.0, 7.0),
    "His148 HD1": (1.0, 6.0),
    "Tyr145 HH": (1.0, 6.0),
    "Thr64 O": (2.5, 6.5),
    "Gln94 NE2": (2.5, 6.5),
}

paths_system = {
    "Reduced": "005-rogfp-glh-md",
    "Oxidized": "007-rogfp-oxd-glh-md",
    "Cu(I)": "006-rogfp-cu-glh-md",
    "Na+": "008-rogfp-na-glh-md",
}
labels_sys_order = ["Reduced", "Oxidized", "Na+", "Cu(I)"]
colors_sys = {
    "Reduced": "#1e2e79",
    "Oxidized": "#EC4067",
    "Cu(I)": "#f99752",
    "Na+": "#1b998b",
}


def populate_data(
    paths_system, paths_data, x_bounds=(0.5, 15.0), bin_width=0.05, bw_method=0.04
):
    data = {}

    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    y_maxes = {}

    for label_data, path_dist in paths_data.items():
        data[label_data] = {}
        y_maxes[label_data] = 0.0
        for label_sys, path_sys in paths_system.items():
            file_path = os.path.join(
                base_dir, f"analysis/{path_sys}/data/struct-desc/{path_dist}.npy"
            )
            dists = np.load(file_path)
            kde = gaussian_kde(dists, bw_method=bw_method)
            y_values = kde(x_values)
            y_max = np.max(y_values)
            if y_max > y_maxes[label_data]:
                y_maxes[label_data] = y_max

            data[label_data][label_sys] = y_values
    return x_values, data, y_maxes


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    n_systems = len(paths_system.keys())
    n_dists = len(paths_data.keys())
    x_values, data, y_maxes = populate_data(paths_system, paths_data)

    gs = grid_spec.GridSpec(nrows=n_systems, ncols=n_dists, bottom=0.14, top=1.09)
    fig = plt.figure(figsize=(6.0, 3.0))

    x_lims = (1, 7)
    axes = []
    i_dist = 0
    for label_data, data_dist in data.items():
        i_sys = 0
        for label_sys in labels_sys_order:
            _data = data[label_data][label_sys]
            axes.append(fig.add_subplot(gs[i_sys, i_dist]))
            axes[-1].plot(x_values, _data, color="#FFFFFF", lw=1.0)
            axes[-1].fill_between(
                x_values, _data, alpha=1, color=colors_sys[label_sys], lw=0.0
            )

            axes[-1].set_xlim(*x_lims_data[label_data])
            rect = axes[-1].patch
            rect.set_alpha(0)
            axes[-1].set_ylim(0, y_maxes[label_data])
            axes[-1].set_yticks([])
            axes[-1].set_yticklabels([])

            spines = ["top", "right", "left", "bottom"]
            for s in spines:
                axes[-1].spines[s].set_visible(False)

            if i_sys % n_systems != n_systems - 1:
                axes[-1].set_xticks([])
                axes[-1].set_xticklabels([])
            else:
                axes[-1].set_xlabel(label_data)

            if i_dist == 0:
                axes[-1].text(
                    x_lims_data[label_data][0] - 0.2,
                    0,
                    label_sys,
                    fontweight="bold",
                    ha="right",
                )

            i_sys += 1
        i_dist += 1

    gs.update(hspace=-0.75)

    fig.savefig("cro66-interaction-ridgeline.svg")
