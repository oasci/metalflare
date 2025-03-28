#!/usr/bin/env python3

import os

import MDAnalysis as mda
from MDAnalysis import transformations

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

    data_dir = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/traj")
    os.makedirs(data_dir, exist_ok=True)
    path_pdb = os.path.join(data_dir, "povme.pdb")

    u = mda.Universe(topology_path, trajectory_paths)
    atoms_str = "protein and (not resname CRO)"
    atoms = u.select_atoms(atoms_str)

    ref_u = u.copy()
    reference = ref_u.select_atoms(atoms_str)

    ag = u.atoms
    workflow = (
        transformations.unwrap(ag),
        transformations.center_in_box(atoms, center="mass"),
        transformations.fit_rot_trans(atoms, reference),
    )
    u.trajectory.add_transformations(*workflow)

    stride = 1
    with mda.Writer(path_pdb, atoms.n_atoms) as W:
        for ts in u.trajectory[None:None:stride]:
            W.write(u.select_atoms(atoms_str))


if __name__ == "__main__":
    main()
