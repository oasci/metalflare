#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rdf import InterRDF


def generate_trajectory_paths(base_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/003-rogfp-cu-md/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
                )
                for prod in range(*prod_range)
            ]
        )
    return [
        path for path in trajectory_paths if "run-01/outputs/10_prod_npt" not in path
    ]


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_paths = generate_trajectory_paths(base_dir)
    topology_path = os.path.join(
        base_dir, "data/003-rogfp-cu-md/simulations/02-prep/mol.prmtop"
    )
    cro_oh_str = "resname CRO and name OH"

    data_dir = os.path.join(base_dir, "analysis/003-rogfp-cu-md/data/rdf")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    cro_oh = u.select_atoms(cro_oh_str)

    for atom_type in ["H", "O", "N"]:
        print(f"Working on {atom_type}")
        atoms = u.select_atoms(f"not resname CRO and element {atom_type}")
        rdf = InterRDF(cro_oh, atoms, nbins=150, range=(0.0, 15.0), norm="rdf")
        rdf.run(step=1)

        data_path = os.path.join(data_dir, f"cro_{atom_type.lower()}_rdf")
        np.save(data_path + "_bins.npy", rdf.results.bins)
        np.save(data_path + "_density.npy", rdf.results.rdf)


if __name__ == "__main__":
    main()
