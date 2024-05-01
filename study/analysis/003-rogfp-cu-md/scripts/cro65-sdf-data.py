#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.density import DensityAnalysis


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_path = os.path.join(
        base_dir, "analysis/003-rogfp-cu-md/data/traj/aligned_cro.nc"
    )
    topology_path = os.path.join(
        base_dir, "data/003-rogfp-cu-md/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, "analysis/003-rogfp-cu-md/data/sdf")
    os.makedirs(data_dir, exist_ok=True)

    element_types = ["O", "N"]
    resids = ["143", "146", "201", "220"]
    for element_type in element_types:
        for resid in resids:
            print(f"Working on resid{resid}_{element_type}")
            u = mda.Universe(topology_path, trajectory_path)
            atoms_1_str = f"(resname CRO and name OH) or (element {element_type} and (resid {resid}))"
            atoms_1 = u.select_atoms(atoms_1_str)
            D = DensityAnalysis(atoms_1, delta=0.5)
            D.run(step=1)

            data_path = os.path.join(
                data_dir, f"resid65_oh-resid{resid}_{element_type.lower()}"
            )
            D.results.density.grid /= D.results.density.grid.sum()
            D.results.density.export(data_path + ".dx")


if __name__ == "__main__":
    main()
