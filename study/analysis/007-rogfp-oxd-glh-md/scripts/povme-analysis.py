#!/usr/bin/env python3

import os
import shutil

from povme import POVME
from povme import enable_logging

enable_logging(20)

SIM_LABEL = os.path.dirname(os.path.abspath(__file__)).split("/")[-2]

def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    path_pdb = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/traj/povme.pdb")
    path_config = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/pocket/povme.yml")

    dir_output = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/pocket/povme/")
    if os.path.exists(dir_output):
        shutil.rmtree(dir_output)

    pocket_id = POVME(path_config)
    pocket_id.run(path_pdb, dir_output)


if __name__ == "__main__":
    main()
