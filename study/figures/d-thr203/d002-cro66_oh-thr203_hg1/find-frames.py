#!/usr/bin/env python3

import json
import os

import numpy as np

os.chdir(os.path.dirname(os.path.realpath(__file__)))


if __name__ == "__main__":
    base_dir = "../../../"

    desc_to_get = {
        "reduced": {
            "path_data": os.path.join(
                base_dir,
                "analysis/005-rogfp-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
            ),
            "values": [1.75, 3.16, 3.97, 4.47, 5.32, 5.83],
        },
        "oxidized": {
            "path_data": os.path.join(
                base_dir,
                "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
            ),
            "values": [1.75, 3.16, 3.97, 4.47, 5.32, 5.83],
        },
        "cu": {
            "path_data": os.path.join(
                base_dir,
                "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_oh-thr201_hg1-dist.npy",
            ),
            "values": [1.75, 3.16, 3.97, 4.47, 5.32, 5.83],
        },
    }

    json_data = {}
    for key, info in desc_to_get.items():
        json_data[key] = []
        data = np.load(info["path_data"])
        values = np.array(info["values"])
        for value in values:
            data_diff = np.abs(data - value)
            arg_best = np.argmin(data_diff)
            json_data[key].append(
                {"value": str(data[arg_best]), "index": str(arg_best)}
            )

    path_out = "relevant-frames.json"
    with open(path_out, "w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=4)
