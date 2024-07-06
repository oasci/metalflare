#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    rogfp_data_path = os.path.join(
        base_dir, "analysis/001-rogfp-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy"
    )
    rogfp_data = np.load(rogfp_data_path)
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/004-rogfp-oxd-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_cu_data_path = os.path.join(
        base_dir,
        "analysis/003-rogfp-cu-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_data_path)

    mean_rogfp_data = np.nanmean(rogfp_data)
    mean_rogfp_oxd_data = np.nanmean(rogfp_cu_data)
    mean_rogfp_cu_data = np.nanmean(rogfp_cu_data)

    kwargs = {
        "kde": True,
        "stat": "density",
        "fill": True,
        "binwidth": 0.1,
        "kde_kws": {"bw_method": 0.04},
    }
    sns.histplot(rogfp_data, label="Reduced", color="#1e2e79", **kwargs)
    sns.histplot(rogfp_oxd_data, label="Oxidized", color="#EC4067", **kwargs)
    sns.histplot(rogfp_cu_data, label="Cu(I)", color="#f99752", **kwargs)

    plt.xlabel("CRO66 OH - THR203 HG1 Distance [Å]")
    plt.xlim(1, 8)
    plt.ylabel("Density")

    plt.legend()
    plt.tight_layout()

    plt.savefig("b008-cro66_oh-thr203_hg1-hist.svg")
