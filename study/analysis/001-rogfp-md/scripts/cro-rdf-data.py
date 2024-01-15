#!/usr/bin/env python3

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "/ihome/jdurrant/amm503/bgfs/oasci/metalflare/study"
    trajectory_paths = []
    for run_i in range(1, 4):
        trajectory_paths.extend([
            os.path.join(base_dir, f"data/001-rogfp-md/simulations/05-prod/run-0{run_i}/outputs/08_prod_npt.nc"),
            os.path.join(base_dir, f"data/001-rogfp-md/simulations/06-prod/run-0{run_i}/outputs/09_prod_npt.nc"),
            os.path.join(base_dir, f"data/001-rogfp-md/simulations/06-prod/run-0{run_i}/outputs/10_prod_npt.nc"),
        ])
    topology_path = os.path.join(base_dir, "data/001-rogfp-md/simulations/02-prep/mol.prmtop")

    # There was an error with run-01/outputs/10_prod_npt
    trajectory_paths = [path for path in trajectory_paths if "run-01/outputs/10_prod_npt" not in path]
    # Specify the atom selections for the two groups
    cro_oh_str = 'resname CRO and name OH'  # Adjust as needed

    data_dir = os.path.join(base_dir, "analysis/001-rogfp-md/data/rdf")
    os.makedirs(data_dir, exist_ok=True)
    # Load the trajectory and calculate distances
    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    cro_oh = u.select_atoms(cro_oh_str)

    for atom_type in ["C", "H", "O", "N", "S"]:
        print(f"Working on {atom_type}")
        atoms = u.select_atoms(f"not resname CRO and name {atom_type}")
        rdf = InterRDF(cro_oh, atoms, nbins=150, range=(0.0, 15.0), norm="density")
        rdf.run(step=10000)
        
        data_path = os.path.join(data_dir, f"cro_{atom_type.lower()}_rdf")
        np.save(data_path + "_bins.npy", rdf.results.bins)
        np.save(data_path + "_density.npy", rdf.results.rdf)
