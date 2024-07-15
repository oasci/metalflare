#!/usr/bin/env python3

import os
import json

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis import transformations
import numpy as np

os.chdir(os.path.dirname(os.path.realpath(__file__)))

def generate_trajectory_paths(base_dir, dir_name, run_range=(1, 4), prod_range=(8, 11)):
    trajectory_paths = []
    for run_i in range(*run_range):
        trajectory_paths.extend(
            [
                os.path.join(
                    base_dir,
                    f"data/{dir_name}/simulations/05-prod/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc",
                )
                for prod in range(*prod_range)
            ]
        )
    return trajectory_paths

FIG_LABEL = "b002"

if __name__ == "__main__":
    base_dir = "../../../"

    path_json = "relevant-frames.json"
    with open(path_json, "r", encoding="utf-8") as f:
        json_data = json.load(f)
    
    # Write PDBs
    dir_info = {
        "reduced": "005-rogfp-glh-md",
        "oxidized": "007-rogfp-oxd-glh-md",
        "cu": "006-rogfp-cu-glh-md"
    }
    pdb_paths = []
    for label, dirname in dir_info.items():

        trajectory_paths = generate_trajectory_paths(base_dir, dir_name=dirname)
        topology_path = os.path.join(
            base_dir, f"data/{dirname}/simulations/02-prep/mol.prmtop"
        )

        u = mda.Universe(topology_path, trajectory_paths)
        atoms_of_interest = u.select_atoms("protein")
        not_atoms_of_interest = u.select_atoms("not protein")
        transforms = [
            transformations.unwrap(atoms_of_interest),
            transformations.center_in_box(atoms_of_interest),
            transformations.wrap(not_atoms_of_interest),
        ]
        u.trajectory.add_transformations(*transforms)

        for info in json_data[label]:
            value_str = f"{float(info['value']):.3f}"
            frame = int(info["index"])
            u.trajectory[frame]
            path_pdb = f"pdbs/{FIG_LABEL}-{label}-{value_str}.pdb"
            pdb_paths.append(path_pdb)
            atoms = u.select_atoms("protein or resname CRO or resname CU1")
            print(f"Writing {path_pdb}")
            atoms.write(path_pdb, bonds=None)
    
    # Align PDBs
    print("Alinging PDBs")
    u_ref = mda.Universe(pdb_paths[0])
    for pdb_path in pdb_paths[1:]:
        u = mda.Universe(pdb_path)
        align.alignto(u, u_ref, select="protein and name CA and (not resid 227)")
        u.atoms.write(pdb_path, bonds=None)


