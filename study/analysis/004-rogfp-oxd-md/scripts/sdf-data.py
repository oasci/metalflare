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

    element_types = ["O", "N", "H"]
    selections = [
        "resname CRO",
        "resid 143",
        "resid 146",
        "resid 201",
        "resid 203",
        "resid 220",
        "resname WAT and (around 5.0 (resname CRO and name OH)) and (around 5.0 (resid 143 and name OH)) and (around 6.0 (resid 203 and name CB))",
    ]
    for element_type in element_types:
        for selection in selections:
            if "resname WAT" in selection:
                res_label = "watpwire"
                if element_type == "N":
                    continue
            else:
                res_label = selection.replace(" ", "").lower()

            print(f"Working on {selection}_{element_type}")
            u = mda.Universe(topology_path, trajectory_path)
            atoms_str = f"element {element_type} and ({selection})"
            atoms = u.select_atoms(atoms_str)
            D = DensityAnalysis(atoms, delta=0.5)
            D.run(step=1)

            data_path = os.path.join(data_dir, f"{res_label}_{element_type.lower()}")
            D.results.density.export(data_path + ".dx")


if __name__ == "__main__":
    main()
