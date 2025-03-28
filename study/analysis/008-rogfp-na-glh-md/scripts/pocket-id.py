#!/usr/bin/env python3

import os
import shutil

from povme import PocketDetector

SIM_LABEL = os.path.dirname(os.path.abspath(__file__)).split("/")[-2]


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    path_pdb = os.path.join(
        base_dir, f"analysis/{SIM_LABEL}/data/traj/sample-povme.pdb"
    )
    path_config = os.path.join(
        base_dir, f"analysis/{SIM_LABEL}/data/pocket/pocket-id.yml"
    )

    dir_output = os.path.join(
        base_dir, f"analysis/{SIM_LABEL}/data/pocket/sample-povme-id/"
    )
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    pocket_id = PocketDetector(path_config)
    pocket_id.run(path_pdb, dir_output)


if __name__ == "__main__":
    main()
