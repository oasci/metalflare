#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    rogfp_data_path = os.path.join(
        base_dir, "analysis/001-rogfp-md/data/struct-desc/cro65_oh-thr201_cg2-dist.npy"
    )
    rogfp_data = np.load(rogfp_data_path)
    rogfp_cu_data_path = os.path.join(
        base_dir,
        "analysis/003-rogfp-cu-md/data/struct-desc/cro65_oh-thr201_cg2-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_data_path)

    kwargs = {"kde": True, "stat": "density", "fill": True}
    sns.histplot(rogfp_data, label="Unbound", color="#1e2e79", **kwargs)
    sns.histplot(rogfp_cu_data, label="Bound", color="#f99752", **kwargs)

    plt.xlabel("CRO66 OH - THR203 CG2 Distance [Å]")
    plt.xlim(right=6.5)
    plt.ylabel("Density")

    plt.legend()
    plt.tight_layout()

    plt.savefig("b003-cro66_oh-thr203_cg2-hist.svg")
