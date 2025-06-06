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


def compute_dihedralal_angle(coords):
    # Given a set of coordinates, compute the dihedralal angle
    # For example, you can use the dihedral function from MDAnalysis
    return mda.lib.distances.calc_dihedrals(coords)


def main():
    base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
    trajectory_paths = generate_trajectory_paths(base_dir)

    topology_path = os.path.join(
        base_dir, f"data/{SIM_LABEL}/simulations/02-prep/mol.prmtop"
    )

    data_dir = os.path.join(base_dir, f"analysis/{SIM_LABEL}/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    atoms = u.select_atoms("(resid 220 and name CB CG CD OE2)")

    atoms_npy_path = os.path.join(data_dir, "glu220_cb_cg_cd_oe2-dihedral.npy")
    atoms_dihedral_array = np.full((n_frames,), np.nan, dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        coords = atoms.positions
        dihedral_angle = mda.lib.distances.calc_dihedrals(*coords)
        atoms_dihedral_array[i] = dihedral_angle

    np.save(atoms_npy_path, atoms_dihedral_array)

    print(atoms_dihedral_array)


if __name__ == "__main__":
    main()
