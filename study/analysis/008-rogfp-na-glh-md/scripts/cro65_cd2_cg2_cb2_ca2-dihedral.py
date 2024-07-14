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
                    f"data/008-rogfp-na-glh-md/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
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
        base_dir, "data/008-rogfp-na-glh-md/simulations/02-prep/mol.prmtop"
    )
    residue_str = "resname CRO and resid 65"
    atom_names = [
        "CD2",
        "CG2",
        "CB2",
        "CA2",
    ]  # Specify the four atoms on the same residue

    data_dir = os.path.join(base_dir, "analysis/008-rogfp-na-glh-md/data/struct-desc/")
    os.makedirs(data_dir, exist_ok=True)

    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    atoms = u.select_atoms(f"{residue_str} and name {' '.join(atom_names)}")

    atoms_npy_path = os.path.join(data_dir, "cro65_cd2_cg2_cb2_ca2-dihedral.npy")
    atoms_dihedral_array = np.full((n_frames,), np.nan, dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        coords = atoms.positions
        dihedral_angle = mda.lib.distances.calc_dihedrals(*coords)
        atoms_dihedral_array[i] = dihedral_angle
    np.save(atoms_npy_path, atoms_dihedral_array)

    print(atoms_dihedral_array)


if __name__ == "__main__":
    main()
