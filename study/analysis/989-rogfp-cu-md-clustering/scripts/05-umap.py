#!/usr/bin/env python3

import os

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap

os.chdir(os.path.dirname(os.path.realpath(__file__)))

sim_features = "../data/sim-features.parquet"


RANDOM_STATE = 3728921

STRIDE = 50


def main():
    df_features = pd.read_parquet(sim_features)
    labels = df_features["system_label"].values[::STRIDE]
    features = df_features.drop(columns=["system_label"]).values[
        ::STRIDE
    ]  # Removed axis=0

    reducer = umap.UMAP(
        n_neighbors=20,
        min_dist=0.1,
        n_components=2,
        random_state=RANDOM_STATE,
        metric="euclidean",
    )

    embedding = reducer.fit_transform(features)

    # Create a colormap
    colors = []
    for label in labels:
        if label == "rogfp":
            colors.append("#1e2e79")
        elif label == "rogfp_cu":
            colors.append("#f99752")

    plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=3, alpha=0.5)
    plt.savefig("umap.svg")
    print(embedding)


if __name__ == "__main__":
    main()
