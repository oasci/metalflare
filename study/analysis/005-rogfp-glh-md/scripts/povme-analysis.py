#!/usr/bin/env python3

import os
import shutil

from povme.pocket.volume import PocketVolume
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

    povme = PocketVolume(path_config)
    povme.run(path_pdb, output_prefix=dir_output, chunk_size=50)


if __name__ == "__main__":
    main()
