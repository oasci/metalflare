#!/usr/bin/env python3
import os
import json
import numpy as np

def angle_difference(angle1, angle2):
    """Calculate the minimum angle difference accounting for periodicity."""
    return np.minimum(np.abs(angle1 - angle2), 2*np.pi - np.abs(angle1 - angle2))

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    base_dir = "../../../"
    desc_to_get = {
        "reduced": {
            "path_data": os.path.join(
                base_dir,
                "analysis/005-rogfp-glh-md/data/struct-desc/cro65_cg2_cb2_ca2_c2-dihedral.npy",
            ),
            "values": [0, np.pi/4, np.pi/2, (3/4)*np.pi],
        },
        "oxidized": {
            "path_data": os.path.join(
                base_dir,
                "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_cg2_cb2_ca2_c2-dihedral.npy",
            ),
            "values": [0, np.pi/4, np.pi/2, (3/4)*np.pi],
        },
        "cu": {
            "path_data": os.path.join(
                base_dir,
                "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_cg2_cb2_ca2_c2-dihedral.npy",
            ),
            "values": [0, np.pi/4, np.pi/2, (3/4)*np.pi],
        },
    }

    json_data = {}
    for key, info in desc_to_get.items():
        json_data[key] = []
        data = np.load(info["path_data"])
        values = np.array(info["values"])

        for value in values:
            # Use the angle_difference function to account for periodicity
            data_diff = angle_difference(data, value)
            arg_best = np.argmin(data_diff)
            json_data[key].append(
                {"value": str(data[arg_best]), "index": str(arg_best)}
            )

    path_out = "relevant-frames.json"
    with open(path_out, "w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=4)