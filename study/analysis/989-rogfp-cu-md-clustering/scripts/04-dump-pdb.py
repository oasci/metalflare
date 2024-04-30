#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align

os.chdir(os.path.dirname(os.path.realpath(__file__)))

representatives_csv_path = "../data/cluster-representatives.csv"


def generate_trajectory_paths(base_dir, sim_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    sim_dir,
                    f"run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
                )
                for prod in range(*prod_range)
            ]
        )
    return trajectory_paths


def main():
    df = pd.read_csv(representatives_csv_path)

    rogfp2_path = "data/001-rogfp-md/simulations"
    rogfp2_cu_path = "data/003-rogfp-cu-md/simulations"

    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    traj_paths_rogfp2 = generate_trajectory_paths(base_dir, f"{rogfp2_path}/05-prod")
    traj_paths_rogfp2_cu = generate_trajectory_paths(
        base_dir, f"{rogfp2_cu_path}/05-prod"
    )

    topology_rogfp2_path = os.path.join(base_dir, f"{rogfp2_path}/02-prep/mol.prmtop")
    topology_rogfp2_cu_path = os.path.join(
        base_dir, f"{rogfp2_cu_path}/02-prep/mol.prmtop"
    )

    universes = {
        "rogfp2": mda.Universe(topology_rogfp2_path, traj_paths_rogfp2),
        "rogfp2_cu": mda.Universe(topology_rogfp2_cu_path, traj_paths_rogfp2_cu),
    }
    pdb_dirs = "../pdbs"
    cluster_i = 0
    atoms_ref = universes["rogfp2"].atoms.select_atoms("backbone")
    for sim_label, idx in df.itertuples(index=False):
        universes[sim_label].trajectory[idx]
        atoms_idx = universes[sim_label].atoms.select_atoms("protein or resname CRO")
        results = align.alignto(atoms_idx.atoms.select_atoms("backbone"), atoms_ref)

        with mda.Writer(os.path.join(pdb_dirs, f"{cluster_i}-{sim_label}.pdb")) as w:
            w.write(atoms_idx)
        cluster_i += 1


if __name__ == "__main__":
    main()
