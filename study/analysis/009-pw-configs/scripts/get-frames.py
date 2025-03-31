#!/usr/bin/env python3

import json
import os

import MDAnalysis as mda
from MDAnalysis import transformations

# === CONFIGURATION === #
base_dir = "/ihome/jdurrant/amm503/ix/oasci/metalflare/study"
cluster_json_path = os.path.join(
    base_dir, "analysis/009-pw-configs/data/cluster-center-closest-frames.json"
)
output_dir = os.path.join(base_dir, "analysis/009-pw-configs/data/cluster_pdbs")
os.makedirs(output_dir, exist_ok=True)

topology_paths = {
    "reduced": os.path.join(
        base_dir, "data/005-rogfp-glh-md/simulations/02-prep/mol.prmtop"
    ),
    "oxidized": os.path.join(
        base_dir, "data/007-rogfp-oxd-glh-md/simulations/02-prep/mol.prmtop"
    ),
    "cu": os.path.join(
        base_dir, "data/006-rogfp-cu-glh-md/simulations/02-prep/mol.prmtop"
    ),
}

trajectory_dirs = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}


def generate_trajectory_paths(
    base_dir, system_name, run_range=(1, 4), prod_range=(8, 11)
):
    sim_base = f"data/{system_name}/simulations/05-prod"
    trajectory_paths = []
    for run_i in range(*run_range):
        for prod in range(*prod_range):
            path = os.path.join(
                base_dir, f"{sim_base}/run-0{run_i}/outputs/{prod:02d}_prod_npt.nc"
            )
            trajectory_paths.append(path)
    return trajectory_paths


def load_cluster_assignments(path):
    with open(path, "r") as f:
        return json.load(f)


def extract_and_write_frames_for_condition(condition, cluster_frames):
    topo = topology_paths[condition]
    system_dir = trajectory_dirs[condition]
    traj_paths = generate_trajectory_paths(base_dir, system_dir)

    print(f"Loading trajectory for {condition} ({len(cluster_frames)} cluster frames)")
    u = mda.Universe(topo, traj_paths)

    sel_str = "protein or resname CRO"
    atoms_of_interest = u.select_atoms(sel_str)
    not_atoms_of_interest = u.select_atoms(f"not ({sel_str})")
    transforms = [
        transformations.unwrap(atoms_of_interest),
        transformations.center_in_box(atoms_of_interest, wrap=True),
        transformations.wrap(not_atoms_of_interest),
    ]
    u.trajectory.add_transformations(*transforms)

    for cluster_label, frame_index in cluster_frames:
        print(f"Extracting frame {frame_index} for cluster {cluster_label}")
        u.trajectory[frame_index]
        out_path = os.path.join(output_dir, f"cluster_{cluster_label}.pdb")
        print(f"Writing {out_path}")
        u.atoms.write(out_path, bonds=None)


def main():
    cluster_assignments = load_cluster_assignments(cluster_json_path)

    # Group clusters by condition
    condition_clusters = {"reduced": [], "oxidized": [], "cu": []}  # type: ignore

    for label, entry in cluster_assignments.items():
        condition_clusters[entry["condition"]].append((label, entry["frame"]))

    for condition, cluster_frames in condition_clusters.items():
        if cluster_frames:
            extract_and_write_frames_for_condition(condition, cluster_frames)


if __name__ == "__main__":
    main()
