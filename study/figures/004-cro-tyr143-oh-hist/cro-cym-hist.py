#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../"

    rogfp_dist_path = os.path.join(
        base_dir, "analysis/001-rogfp-md/data/cro_thr143_oh.npy"
    )
    rogfp_dist = np.load(rogfp_dist_path)
    rogfp_cu_dist_path = os.path.join(
        base_dir, "analysis/003-rogfp-cu-md/data/cro_thr143_oh.npy"
    )
    rogfp_cu_dist = np.load(rogfp_cu_dist_path)

    mean_rogfp_dist = np.nanmean(rogfp_dist)
    mean_rogfp_cu_dist = np.nanmean(rogfp_cu_dist)

    print(f"Mean rogfp dist:    {mean_rogfp_dist:.3f}")  # 7.463
    print(f"Mean rogfp cu dist: {mean_rogfp_cu_dist:.3f}")  # 6.605

    kwargs = {"kde": True, "stat": "density", "fill": True}
    sns.histplot(rogfp_dist, label="Unbound", color="#1e2e79", **kwargs)
    sns.histplot(rogfp_cu_dist, label="Bound", color="#f99752", **kwargs)

    plt.xlabel("OH CRO65 - OH THR143 Distance [Å]")
    plt.ylabel("Density")

    plt.legend()

    plt.savefig("cro-thr-hist.png")
