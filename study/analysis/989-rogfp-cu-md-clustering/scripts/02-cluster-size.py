#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

os.chdir(os.path.dirname(os.path.realpath(__file__)))

sim_features = "../data/sim-scaled-features.parquet"


RANDOM_STATE = 3728921

STRIDE = 50


def main():
    df_features = pd.read_parquet(sim_features)
    data_all = df_features.drop(columns=["system_label"], axis=0).values[::STRIDE]

    results = {"n_clusters": [], "silhouette_score": []}
    n_clusters_range = tuple(range(4, 30))
    for n_clusters in n_clusters_range:
        kmeans_model = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE).fit(
            data_all
        )
        silhouette_avg = silhouette_score(data_all, kmeans_model.labels_)
        results["n_clusters"].append(n_clusters)
        results["silhouette_score"].append(silhouette_avg)
        print(
            f"n_clusters = {n_clusters} average silhouette_score is : {silhouette_avg}"
        )

    print(results)


if __name__ == "__main__":
    main()
