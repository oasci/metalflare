#!/usr/bin/env python3

import os

import MDAnalysis as mda
import numpy as np
import pandas as pd

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

    base_dir = "/ihome/jdurrant/amm503/bgfs/oasci/metalflare/study"
    traj_paths_rogfp2 = generate_trajectory_paths(base_dir, f"{rogfp2_path}/05-prod")
    traj_paths_rogfp2_cu = generate_trajectory_paths(
        base_dir, f"{rogfp2_cu_path}/05-prod"
    )

    topology_rogfp2_path = os.path.join(base_dir, f"{rogfp2_path}/02-prep/mol.prmtop")
    topology_rogfp2_cu_path = os.path.join(
        base_dir, f"{rogfp2_cu_path}/02-prep/mol.prmtop"
    )

    u_rogfp2 = mda.Universe(topology_rogfp2_path, traj_paths_rogfp2)
    u_rogfp2_cu = mda.Universe(topology_rogfp2_cu_path, traj_paths_rogfp2_cu)


if __name__ == "__main__":
    main()
