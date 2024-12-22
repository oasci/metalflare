#!/usr/bin/env python3

import os

import numpy as np

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

    red_data_path = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65-rmsf.npy",
    )
    red_data = np.load(red_data_path)

    # Oxidized
    oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65-rmsf.npy",
    )
    oxd_data = np.load(oxd_data_path)

    cu_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65-rmsf.npy",
    )
    cu_data = np.load(cu_path)

    na_path = os.path.join(
        base_dir,
        "analysis/008-rogfp-na-glh-md/data/struct-desc/cro65-rmsf.npy",
    )
    na_data = np.load(na_path)
