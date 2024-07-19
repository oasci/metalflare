#!/usr/bin/env python3

import os

import MDAnalysis as mda
from MDAnalysis import transformations
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
    atoms1_str = "resid 203 and name HG"
    atoms2_str = "resname WAT and (around 6.0 (resname CRO and name OH)) and (around 6.5 (resid 143 and name OH)) and (around 6.0 (resid 201 and name CA)) and name O"

    data_dir = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    atoms_of_interest = u.select_atoms("protein")
    not_atoms_of_interest = u.select_atoms("not protein")
    transforms = [
        transformations.unwrap(atoms_of_interest),
        transformations.center_in_box(atoms_of_interest, wrap=True),
        transformations.wrap(not_atoms_of_interest),
    ]
    u.trajectory.add_transformations(*transforms)
    n_frames = len(u.trajectory)

    atoms_npy_path = os.path.join(data_dir, "ser203_hg-h2o_o-dist.npy")
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
