#!/usr/bin/env python3

import os

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
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
    atoms_str = "resname CRO"

    data_dir = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    average = align.AverageStructure(
        u,
        u,
        select=atoms_str,
        ref_frame=0
    ).run()
    ref = average.results.universe
    aligner = align.AlignTraj(
        u, ref,
        select=atoms_str,
    ).run()

    atoms = u.select_atoms(atoms_str)
    R = rms.RMSF(atoms).run()

    atoms_npy_path = os.path.join(data_dir, "cro65-rmsf.npy")
    np.save(atoms_npy_path, R)

    print(R)
    print(f"Min:  {np.min(R)}")
    print(f"Mean: {np.mean(R)}")
    print(f"Max:  {np.max(R)}")


if __name__ == "__main__":
    main()
