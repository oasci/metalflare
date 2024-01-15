#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../"

    for atom_type in ["C", "H", "O", "N", "S"]:
        bins_path = os.path.join(
            base_dir,
            "analysis/001-rogfp-md/data/rdf/",
            f"cro_{atom_type.lower()}_rdf_bins.npy"
        )
        bins = np.load(bins_path)
        density_path = os.path.join(
            base_dir,
            "analysis/001-rogfp-md/data/rdf/",
            f"cro_{atom_type.lower()}_rdf_density.npy"
        )
        density = np.load(density_path)
        
        kwargs = {"kde": True, "stat": "density", "fill": True, "bins": bins}
        sns.histplot(density, **kwargs)

        plt.xlabel("Distance [Å]")
        plt.ylabel(f"g(O, {atom_tpye})")

        plt.savefig("cro-cym-hist.png")
