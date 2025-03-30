#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pes import create_pes

os.chdir(os.path.dirname(os.path.realpath(__file__)))

fig_label = "fig003"
data_x_str = "cys202_c-ser203_n_ca_c-dihedral"
data_x_label = "Ser205 $\phi$ [°]"
data_y_str = "ser203_n_ca_c-ala204_n-dihedral"
data_y_label = "Ser205 $\psi$ [°]"

bin_width_angle = 2
bin_width_dist = 0.1
bins = (
    np.arange(-360, 0 + 0.01, bin_width_angle),
    np.arange(0, 360 + 0.01, bin_width_angle),
)
x_lims = (-210, -30)
x_ticks = np.arange(-210, 0 + 0.001, 60)
y_lims = (60, 220)
y_ticks = np.arange(min(y_lims), max(y_lims) + 0.001, 30)

pes_vmin = 0
pes_vmax = 3.7
levels = 8
T = 300.0

label_x = 0.5
label_y = 0.965

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
use_mpl_rc_params(rc_json_path, font_dirs)

# Load and process original data for each case

# Reduced data
path_x_red = os.path.join(
    base_dir,
    f"analysis/005-rogfp-glh-md/data/struct-desc/{data_x_str}.npy",
)
data_x_red = np.load(path_x_red)
data_x_red = np.degrees(data_x_red)
data_x_red[data_x_red >= 0] -= 360

path_y_red = os.path.join(
    base_dir,
    f"analysis/005-rogfp-glh-md/data/struct-desc/{data_y_str}.npy",
)
data_y_red = np.load(path_y_red)
data_y_red = np.degrees(data_y_red)
data_y_red[data_y_red <= 0] += 360

# Oxidized data
path_x_oxd = os.path.join(
    base_dir,
    f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_x_str}.npy",
)
data_x_oxd = np.load(path_x_oxd)
data_x_oxd = np.degrees(data_x_oxd)
data_x_oxd[data_x_oxd >= 0] -= 360

path_y_oxd = os.path.join(
    base_dir,
    f"analysis/007-rogfp-oxd-glh-md/data/struct-desc/{data_y_str}.npy",
)
data_y_oxd = np.load(path_y_oxd)
data_y_oxd = np.degrees(data_y_oxd)
data_y_oxd[data_y_oxd <= 0] += 360

# Cu data
path_x_cu = os.path.join(
    base_dir,
    f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_x_str}.npy",
)
data_x_cu = np.load(path_x_cu)
data_x_cu = np.degrees(data_x_cu)
data_x_cu[data_x_cu >= 0] -= 360

path_y_cu = os.path.join(
    base_dir,
    f"analysis/006-rogfp-cu-glh-md/data/struct-desc/{data_y_str}.npy",
)
data_y_cu = np.load(path_y_cu)
data_y_cu = np.degrees(data_y_cu)
data_y_cu[data_y_cu <= 0] += 360

# Create a single figure with three horizontal subplots
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.5))

# Plot for the reduced dataset
axis_i = 0
create_pes(
    data_x_red,
    data_y_red,
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
for y_tick in y_ticks[1:]:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for x_tick in x_ticks[1:-1]:
    axes[axis_i].axvline(x_tick, **kwargs_grid)


# Plot for the oxidized dataset
axis_i = 1
create_pes(
    data_x_oxd,
    data_y_oxd,
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
for y_tick in y_ticks[1:]:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for x_tick in x_ticks[1:-1]:
    axes[axis_i].axvline(x_tick, **kwargs_grid)


# Plot for the Cu dataset
axis_i = 2
create_pes(
    data_x_cu,
    data_y_cu,
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
for y_tick in y_ticks[1:]:
    axes[axis_i].axhline(y_tick, **kwargs_grid)
for x_tick in x_ticks[1:-1]:
    axes[axis_i].axvline(x_tick, **kwargs_grid)


# Adding color bar

divider = make_axes_locatable(axes[axis_i])
cax = divider.append_axes("right", size="4%", pad=0.00)
fig.colorbar(
    axes[axis_i].collections[0],  # use the contour set from the Cu plot
    cax=cax,
    label="PMF [kcal/mol]",
    ticks=list(range(int(pes_vmin), int(np.floor(pes_vmax)) + 1)),
)
axes[axis_i].axhline(
    180,
    linestyle=(0, (0.1, 3)),
    linewidth=1.5,
    color=color_grid,
    zorder=0,
    dash_capstyle="round",
)
axes[axis_i].axhline(
    120,
    linestyle=(0, (0.1, 3)),
    linewidth=1.5,
    color=color_grid,
    zorder=0,
    dash_capstyle="round",
)

plt.tight_layout()
fig.savefig(f"{fig_label}.svg")
plt.close(fig)
