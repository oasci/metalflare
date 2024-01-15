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
        density_path = os.path.join(
            base_dir,
            "analysis/001-rogfp-md/data/rdf/",
            f"cro_{atom_type.lower()}_rdf_density.npy",
        )
        density = np.load(density_path)

        kwargs = {
            "kde": True,
            "stat": "density",
            "fill": True,
            "binrange": (0, 15),
            "bins": 150,
        }
        sns.histplot(density, **kwargs)

        plt.xlabel("Distance [Å]")
        plt.ylabel(f"g(O, {atom_type})")

        plt.savefig(f"cro-rdf-{atom_type}.png")
        plt.clf()
