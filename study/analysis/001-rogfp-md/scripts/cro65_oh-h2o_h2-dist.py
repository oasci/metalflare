#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np


def generate_trajectory_paths(base_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/001-rogfp-md/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
                )
                for prod in range(*prod_range)
            ]
        )
    return trajectory_paths


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_paths = generate_trajectory_paths(base_dir)

    topology_path = os.path.join(
        base_dir, "data/001-rogfp-md/simulations/02-prep/mol.prmtop"
    )
    atoms1_str = "resname CRO and name OH"
    atoms2_str = "resname WAT and (around 6.0 (resname CRO and name OH)) and (around 6.5 (resid 143 and name OH)) and (around 6.0 (resid 201 and name CA)) and name H2"

    data_dir = os.path.join(base_dir, "analysis/001-rogfp-md/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    atoms_npy_path = os.path.join(data_dir, "cro65_oh-h2o_h2-dist.npy")
    atoms_dist_array = np.full((n_frames,), np.nan, dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        atoms_1 = u.select_atoms(atoms1_str)
        atoms_2 = u.select_atoms(atoms2_str)
        try:
            dist = np.linalg.norm(atoms_1.positions - atoms_2.positions)
        except:
            dist = np.nan  # Dummy value when no water present
        atoms_dist_array[i] = dist

    np.save(atoms_npy_path, atoms_dist_array)

    print(atoms_dist_array)


if __name__ == "__main__":
    main()
