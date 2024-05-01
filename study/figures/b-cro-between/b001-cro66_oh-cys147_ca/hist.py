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
        base_dir, "analysis/001-rogfp-md/data/struct-desc/cro65_oh-cym145_ca-dist.npy"
    )
    rogfp_dist = np.load(rogfp_dist_path)
    rogfp_cu_dist_path = os.path.join(
        base_dir,
        "analysis/003-rogfp-cu-md/data/struct-desc/cro65_oh-cym145_ca-dist.npy",
    )
    rogfp_cu_dist = np.load(rogfp_cu_dist_path)

    mean_rogfp_dist = np.nanmean(rogfp_dist)
    mean_rogfp_cu_dist = np.nanmean(rogfp_cu_dist)

    bw_method = 0.1
    kwargs = {
        "kde": True,
        "stat": "density",
        "fill": True,
        "kde_kws": {"bw_method": bw_method},
    }
    sns.histplot(rogfp_dist, label="Unbound", color="#1e2e79", **kwargs)
    sns.histplot(rogfp_cu_dist, label="Bound", color="#f99752", **kwargs)

    plt.xlabel(r"CRO66 OH - CYS147 CA Distance [Å]")
    plt.xlim(4.5, 9)
    plt.ylabel("Density")

    plt.legend()
    plt.tight_layout()

    plt.savefig("b001-cro66_oh-cys147_ca-hist.svg")
