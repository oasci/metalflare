#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align

os.chdir(os.path.dirname(os.path.realpath(__file__)))

representatives_csv_path = "../data/representative_structures.csv"


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
        "rogfp": mda.Universe(topology_rogfp2_path, traj_paths_rogfp2),
        "rogfp_cu": mda.Universe(topology_rogfp2_cu_path, traj_paths_rogfp2_cu),
    }

    # Save all unqiue structures
    single_select_df = df[df.isna().sum(axis=1) == 6]
    pdb_dirs = "../pdbs/single-select"
    os.makedirs(pdb_dirs, exist_ok=True)
    single_select_df.to_csv(os.path.join(pdb_dirs, "single-select.csv"), index=False)
    structure_i = 0
    for _, row in single_select_df.iterrows():
        sim_label = row["system_label"]
        idx = row["traj_idx"]
        universes[sim_label].trajectory[idx]
        atoms_idx = universes[sim_label].atoms.select_atoms(
            "protein or resname CRO or resname CU1"
        )

        with mda.Writer(
            os.path.join(pdb_dirs, f"single-select-{structure_i}-{sim_label}.pdb")
        ) as w:
            w.write(atoms_idx)
        structure_i += 1


if __name__ == "__main__":
    main()
