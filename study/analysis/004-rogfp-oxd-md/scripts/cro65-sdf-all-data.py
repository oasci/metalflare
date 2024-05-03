#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.density import DensityAnalysis


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_path = os.path.join(
        base_dir, "analysis/004-rogfp-oxd-md/data/traj/aligned_cro.nc"
    )
    topology_path = os.path.join(
        base_dir, "data/004-rogfp-oxd-md/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, "analysis/004-rogfp-oxd-md/data/sdf")
    os.makedirs(data_dir, exist_ok=True)

    element_types = ["O", "N"]
    resids = ["143", "146", "201", "203", "220"]
    for element_type in element_types:
        print(f"Working on element {element_type}")
        resid_string = "resid " + " or resid ".join([resid for resid in resids])
        u = mda.Universe(topology_path, trajectory_path)
        atoms_1_str = f"(resname CRO and name OH) or (element {element_type} and ({resid_string}))"
        if element_type == "O":
            atoms_1_str += "or (resname WAT and (around 5.0 (resname CRO and name OH)) and (around 5.0 (resid 143 and name OH)) and (around 6.0 (resid 203 and name CB)))"
        atoms_1 = u.select_atoms(atoms_1_str)
        D = DensityAnalysis(atoms_1, delta=0.5)
        D.run(step=1)

        D.results.density.grid /= D.results.density.grid.sum()

        data_path = os.path.join(data_dir, f"resid65_oh-{element_type.lower()}")
        D.results.density.export(data_path + ".dx")


if __name__ == "__main__":
    main()
