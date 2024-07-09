#!/usr/bin/env python3

import os
from random import randrange

import MDAnalysis as mda
from MDAnalysis import transformations


def generate_trajectory_paths(base_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/005-rogfp-glh-md/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
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
        base_dir, "data/005-rogfp-glh-md/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, "analysis/005-rogfp-glh-md/data/traj")
    os.makedirs(data_dir, exist_ok=True)

    for _ in range(5):
        u = mda.Universe(topology_path, trajectory_paths)
        atoms_of_interest = u.select_atoms("protein")
        not_atoms_of_interest = u.select_atoms("not protein")

        n_max = len(u.trajectory)

        i_selection = randrange(n_max)
        path_pdb = os.path.join(data_dir, f"frame_{i_selection}.pdb")

        u.trajectory[i_selection]
        transforms = [
            transformations.unwrap(atoms_of_interest),
            transformations.center_in_box(atoms_of_interest, wrap=True),
            transformations.wrap(not_atoms_of_interest),
        ]
        u.trajectory.add_transformations(*transforms)
        print(f"Writing {path_pdb}")
        u.atoms.write(path_pdb, bonds=None)


if __name__ == "__main__":
    main()
