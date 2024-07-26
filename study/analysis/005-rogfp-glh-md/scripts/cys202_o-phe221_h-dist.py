#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np

SIM_LABEL = os.path.dirname(os.path.abspath(__file__)).split("/")[-2]


def generate_trajectory_paths(base_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/{SIM_LABEL}/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
                )
                for prod in range(*prod_range)
            ]
        )
    return trajectory_paths


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_paths = generate_trajectory_paths(base_dir)

    topology_path = os.path.join(
        base_dir, f"data/{SIM_LABEL}/simulations/02-prep/mol.prmtop"
    )
    atoms1_str = "resid 202 and name O"
    atoms2_str = "resid 221 and name H"

    data_dir = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    atoms_1 = u.select_atoms(atoms1_str)
    atoms_2 = u.select_atoms(atoms2_str)

    atoms_npy_path = os.path.join(data_dir, "cys202_o-phe221_h-dist.npy")
    atoms_dist_array = np.full((n_frames,), np.nan, dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        dist = np.linalg.norm(atoms_1.positions - atoms_2.positions)
        atoms_dist_array[i] = dist
    np.save(atoms_npy_path, atoms_dist_array)

    print(atoms_dist_array)


if __name__ == "__main__":
    main()
