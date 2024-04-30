#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap

os.chdir(os.path.dirname(os.path.realpath(__file__)))

sim_features = "../data/sim-features.parquet"


STRIDE = 1
USE_TMP = False


def main():
    df_features = pd.read_parquet(sim_features)
    labels = df_features["system_label"].values[::STRIDE]
    features = df_features.drop(columns=["system_label"]).values[::STRIDE]

    kwargs = {
        "n_neighbors": 15,
        "n_components": 2,
        "min_dist": 0.0,
        "densmap": False,
        "dens_lambda": 2.0,
    }

    reducer = umap.UMAP(**kwargs)

    embedding = reducer.fit_transform(features)

    # Create a colormap
    colors = []
    for label in labels:
        if label == "rogfp_red":
            colors.append("#1e2e79")
        elif label == "rogfp_oxd":
            colors.append("#ec4067")
        elif label == "rogfp_cu":
            colors.append("#f99752")

    plt.scatter(
        embedding[:, 0], embedding[:, 1], c=colors, s=50, alpha=0.4, edgecolors="none"
    )

    if USE_TMP:
        save_dir = "../tmp"
    else:
        save_dir = "../data"
    umap_name = f"umap-neigh{kwargs['n_neighbors']}-dist{kwargs['min_dist']}"
    if reducer.densmap:
        umap_name += "-dens"
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{umap_name}.png", dpi=1000)
    plt.close()

    np.save(f"{save_dir}/{umap_name}-embeddings.npy", embedding)


if __name__ == "__main__":
    main()
