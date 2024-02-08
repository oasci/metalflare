#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

os.chdir(os.path.dirname(os.path.realpath(__file__)))

rogfp2_file_path = "../data/rogfp2-features.npy"
rogfp2_cu_file_path = "../data/rogfp2-cu-features.npy"


RANDOM_STATE = 3728921

STRIDE = 10


def main():
    rogfp2_arr = np.load(rogfp2_file_path)
    rogfp2_cu_arr = np.load(rogfp2_cu_file_path)

    rogfp2_arr = rogfp2_arr[::STRIDE]
    rogfp2_cu_arr = rogfp2_cu_arr[::STRIDE]

    simulation_labels = np.empty(
        rogfp2_arr.shape[0] + rogfp2_cu_arr.shape[0], dtype=np.uint8
    )
    simulation_labels[: rogfp2_arr.shape[0]] = 0
    simulation_labels[rogfp2_arr.shape[0] :] = 1

    data_all = np.concatenate((rogfp2_arr, rogfp2_cu_arr))

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
    # n_clusters = 7 average silhouette_score is : 0.32744956151959825


if __name__ == "__main__":
    main()
