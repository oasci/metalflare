#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_path = os.path.join(
        base_dir, "analysis/004-rogfp-oxd-md/data/traj/aligned_cro.nc"
    )
    topology_path = os.path.join(
        base_dir, "data/004-rogfp-oxd-md/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, "analysis/004-rogfp-oxd-md/data/structure")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_path)
    data_path = os.path.join(data_dir, "average-structure.pdb")
    align.AverageStructure(
        u, filename=data_path, select="protein or resname CRO", ref_frame=0
    ).run()


if __name__ == "__main__":
    main()
