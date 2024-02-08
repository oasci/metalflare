#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

os.chdir(os.path.dirname(os.path.realpath(__file__)))

rogfp2_file_path = "../data/rogfp2-features.npy"
rogfp2_cu_file_path = "../data/rogfp2-cu-features.npy"


RANDOM_STATE = 3728921
STRIDE = 1
N_CLUSTERS = 7


def main():
    rogfp2_arr = np.load(rogfp2_file_path)
    rogfp2_cu_arr = np.load(rogfp2_cu_file_path)

    rogfp2_arr = rogfp2_arr[::STRIDE]
    rogfp2_cu_arr = rogfp2_cu_arr[::STRIDE]

    idx_split = rogfp2_arr.shape[0]

    simulation_labels = np.empty(
        rogfp2_arr.shape[0] + rogfp2_cu_arr.shape[0], dtype=np.uint8
    )
    simulation_labels[:idx_split] = 0
    simulation_labels[idx_split:] = 1

    data_all = np.concatenate((rogfp2_arr, rogfp2_cu_arr))

    kmeans_model = KMeans(n_clusters=N_CLUSTERS, random_state=RANDOM_STATE).fit(
        data_all
    )
    cluster_labels = kmeans_model.labels_
    closest, _ = pairwise_distances_argmin_min(kmeans_model.cluster_centers_, data_all)

    np.save("../data/cluster-labels-rogfp2.npy", cluster_labels[simulation_labels == 0])
    np.save(
        "../data/cluster-labels-rogfp2-cu.npy", cluster_labels[simulation_labels == 1]
    )

    rep_structures = {"simulation": [], "step": []}
    for i in closest:
        if i < idx_split:
            sim = "rogfp2"
            idx = i
        else:
            sim = "rogfp2_cu"
            idx = i - idx_split
        rep_structures["simulation"].append(sim)
        rep_structures["step"].append(idx)

    pd.DataFrame(rep_structures).to_csv(
        "../data/cluster-representatives.csv", index=False
    )


if __name__ == "__main__":
    main()
