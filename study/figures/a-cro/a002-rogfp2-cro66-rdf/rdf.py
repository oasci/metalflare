#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    for atom_type in ["H", "O", "N"]:
        # bins
        bins_path_unbound = os.path.join(
            base_dir,
            "analysis/001-rogfp-md/data/rdf/",
            f"cro65_{atom_type.lower()}_rdf_bins.npy",
        )
        bin_edges_unbound = np.load(bins_path_unbound)

        bins_path_bound = os.path.join(
            base_dir,
            "analysis/003-rogfp-cu-md/data/rdf/",
            f"cro65_{atom_type.lower()}_rdf_bins.npy",
        )
        bin_edges_bound = np.load(bins_path_bound)

        # rdf
        density_path_unbound = os.path.join(
            base_dir,
            "analysis/001-rogfp-md/data/rdf/",
            f"cro65_{atom_type.lower()}_rdf_density.npy",
        )
        density_unbound = np.load(density_path_unbound)
        density_path_bound = os.path.join(
            base_dir,
            "analysis/003-rogfp-cu-md/data/rdf/",
            f"cro65_{atom_type.lower()}_rdf_density.npy",
        )
        density_bound = np.load(density_path_bound)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(
            bin_edges_unbound,
            density_unbound,
            label="Unbound",
            color="#1e2e79",
            linewidth=2.0,
            zorder=2,
        )
        ax.plot(
            bin_edges_bound,
            density_bound,
            label="Bound",
            color="#f99752",
            linewidth=3.5,
            zorder=1,
        )

        ax.set_xlabel("Distance [Å]")
        ax.set_xlim(0, 10)
        ax.set_ylabel("CRO66 OH $g_{O" + atom_type + "}$")
        ax.set_ylim(bottom=0)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"cro66-rdf-{atom_type.lower()}.svg")
        plt.clf()
