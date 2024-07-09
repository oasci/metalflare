#!/usr/bin/env python3

import os
from random import randrange

import MDAnalysis as mda
from MDAnalysis.transformations import center_in_box, unwrap


def generate_trajectory_paths(base_dir, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/004-rogfp-oxd-md/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
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
        base_dir, "data/004-rogfp-oxd-md/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, "analysis/004-rogfp-oxd-md/data/traj")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)

    n_max = len(u.trajectory)

    i_selection = randrange(n_max)
    path_pdb = os.path.join(data_dir, f"frame_{i_selection}.pdb")
    u.trajectory[i_selection]
    u.trajectory.add_transformations(center_in_box(u.atoms, center="geometry"))
    print(f"Writing {path_pdb}")
    u.atoms.write(path_pdb, bonds=None)


if __name__ == "__main__":
    main()
