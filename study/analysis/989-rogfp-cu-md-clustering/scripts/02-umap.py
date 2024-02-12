#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd
import umap

os.chdir(os.path.dirname(os.path.realpath(__file__)))

sim_features = "../data/sim-features.parquet"


RANDOM_STATE = 3728921

STRIDE = 50


def main():
    df_features = pd.read_parquet(sim_features)
    labels = df_features["system_label"].values[::STRIDE]
    features = df_features.drop(columns=["system_label"]).values[::STRIDE]

    kwargs = {
        "n_neighbors": 20,
        "n_components": 2,
        "min_dist": 0.0,
        "random_state": RANDOM_STATE,
        "densmap": False,
        "dens_lambda": 2.0,
    }

    reducer = umap.UMAP(**kwargs)

    embedding = reducer.fit_transform(features)

    # Create a colormap
    colors = []
    for label in labels:
        if label == "rogfp":
            colors.append("#1e2e79")
        elif label == "rogfp_cu":
            colors.append("#f99752")

    plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=3, alpha=0.5)
    plt.savefig("../umap-tmp/umap.svg")
    plt.close()


if __name__ == "__main__":
    main()
