#!/usr/bin/env python3
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "fig005"
data_x_label = "VAMP 1"
data_y_label = "VAMP 2"

bin_width_x = 0.005
bin_width_y = 0.05
bins = (
    np.arange(-1.3, -0.7 + bin_width_x, bin_width_x),
    np.arange(-2, 4 + bin_width_y, bin_width_y),
)
x_lims = (-1.15, -0.75)
x_ticks = np.arange(-1.2, max(x_lims) + 0.001, 0.1)
y_lims = (-1.6, 4.0)
y_ticks = np.arange(-2, max(y_lims) + 0.001, 1.0)

pes_vmin = 0
pes_vmax = 4.0
levels = 10
T = 300.0

label_x = 0.5
label_y = 0.96

color_grid = "#ededed"
kwargs_grid = {
    "linestyle": (0, (0.1, 3)),
    "linewidth": 1.5,
    "color": color_grid,
    "zorder": 0,
    "dash_capstyle": "round",
}


# Specify the paths to the trajectory and topology files
base_dir = "../../../"

# Update plot params
rc_json_path = os.path.join(base_dir, "misc/003-figure-style/matplotlib-rc-params.json")
font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
use_mpl_rc_params(rc_json_path, font_dirs)  # type: ignore

# Load and process original data for each case

# Reduced data
path_red = os.path.join(
    base_dir,
    "analysis/009-pw-configs/data/reduced-vamp.npy",
)
data_red = np.load(path_red)

# Oxidized data
path_oxd = os.path.join(
    base_dir,
    "analysis/009-pw-configs/data/oxidized-vamp.npy",
)
data_oxd = np.load(path_oxd)

# Cu data
path_cu = os.path.join(
    base_dir,
    "analysis/009-pw-configs/data/cu-vamp.npy",
)
data_cu = np.load(path_cu)

# Load cluster centers with labels (A, B, C, ...)
cluster_json_path = os.path.join(
    base_dir,
    "analysis/009-pw-configs/data/vamp-clustering-global.json",
)
with open(cluster_json_path, "r") as f:
    cluster_centers = json.load(f)["cluster_centers"]

# Convert to array for easier use
cluster_coords = {label: np.array(coords) for label, coords in cluster_centers.items()}


# Create a single figure with three horizontal subplots
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.5))

# Plot for the reduced dataset
axis_i = 0
create_pes(
    data_red[:, 0],
    data_red[:, 1],
    bins=bins,
    vmin=pes_vmin,
    vmax=pes_vmax,
    levels=levels,
    T=T,
    ax=axes[axis_i],
    colorbar=False,
)
axes[axis_i].set_xlabel(data_x_label)
axes[axis_i].set_xlim(*x_lims)
axes[axis_i].set_xticks(x_ticks)
axes[axis_i].set_ylabel(data_y_label)
axes[axis_i].set_ylim(*y_lims)
axes[axis_i].set_yticks(y_ticks)
axes[axis_i].text(
    label_x,
    label_y,
    "Reduced",
    ha="center",
    va="center",
    transform=axes[axis_i].transAxes,
)
for x_tick in x_ticks:
    axes[axis_i].axvline(x_tick, **kwargs_grid)
for y_tick in y_ticks:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for label, coord in cluster_coords.items():
    if x_lims[0] <= coord[0] <= x_lims[1] and y_lims[0] <= coord[1] <= y_lims[1]:
        axes[axis_i].text(
            coord[0],
            coord[1],
            label,
            fontsize=9,
            weight="bold",
            ha="center",
            va="center",
            color="white",
            zorder=10,
        )


# Plot for the oxidized dataset
axis_i = 1
create_pes(
    data_oxd[:, 0],
    data_oxd[:, 1],
    bins=bins,
    vmin=pes_vmin,
    vmax=pes_vmax,
    levels=levels,
    T=T,
    ax=axes[axis_i],
    colorbar=False,
)
axes[axis_i].set_xlabel(data_x_label)
axes[axis_i].set_xlim(*x_lims)
axes[axis_i].set_xticks(x_ticks)
axes[axis_i].set_ylim(*y_lims)
axes[axis_i].set_yticks([])
axes[axis_i].text(
    label_x,
    label_y,
    "Oxidized",
    ha="center",
    va="center",
    transform=axes[axis_i].transAxes,
)
for x_tick in x_ticks:
    axes[axis_i].axvline(x_tick, **kwargs_grid)
for y_tick in y_ticks:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for label, coord in cluster_coords.items():
    if x_lims[0] <= coord[0] <= x_lims[1] and y_lims[0] <= coord[1] <= y_lims[1]:
        axes[axis_i].text(
            coord[0],
            coord[1],
            label,
            fontsize=9,
            weight="bold",
            ha="center",
            va="center",
            color="white",
            zorder=10,
        )


# Plot for the Cu dataset
axis_i = 2
create_pes(
    data_cu[:, 0],
    data_cu[:, 1],
    bins=bins,
    vmin=pes_vmin,
    vmax=pes_vmax,
    levels=levels,
    T=T,
    ax=axes[axis_i],
    colorbar=False,
)
axes[axis_i].set_xlabel(data_x_label)
axes[axis_i].set_xlim(*x_lims)
axes[axis_i].set_xticks(x_ticks)
axes[axis_i].set_ylim(*y_lims)
axes[axis_i].set_yticks([])
axes[axis_i].text(
    label_x,
    label_y,
    "Cu(I)",
    ha="center",
    va="center",
    transform=axes[axis_i].transAxes,
)
for x_tick in x_ticks:
    axes[axis_i].axvline(x_tick, **kwargs_grid)
for y_tick in y_ticks:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for label, coord in cluster_coords.items():
    if x_lims[0] <= coord[0] <= x_lims[1] and y_lims[0] <= coord[1] <= y_lims[1]:
        axes[axis_i].text(
            coord[0],
            coord[1],
            label,
            fontsize=9,
            weight="bold",
            ha="center",
            va="center",
            color="white",
            zorder=10,
        )


# Adding color bar

divider = make_axes_locatable(axes[axis_i])
cax = divider.append_axes("right", size="4%", pad=0.00)
fig.colorbar(
    axes[axis_i].collections[0],  # use the contour set from the Cu plot
    cax=cax,
    label="PMF [kcal/mol]",
    ticks=list(range(int(pes_vmin), int(np.floor(pes_vmax)) + 1)),
)

plt.tight_layout()
fig.savefig(f"{fig_label}.svg")
plt.close(fig)
