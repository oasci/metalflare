#!/usr/bin/env python3

import os

import numpy as np
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

os.chdir(os.path.dirname(os.path.realpath(__file__)))

sim_features = "../data/umap-neigh15-dist0.0-embeddings.npy"


RANDOM_STATE = 3728921

STRIDE = 50


def to_yaml(data, filename):
    """
    Dumps a dictionary `data` into a YAML file named `filename`.

    Parameters:
    - data (dict): The dictionary to be dumped into a YAML file.
    - filename (str): The name of the file to which the dictionary will be dumped.
    """
    with open(filename, "w", encoding="utf-8") as file:
        yaml.dump(data, file, allow_unicode=True)


def main():
    embeddings = np.load(sim_features)[::STRIDE]

    results = {"n_clusters": [], "silhouette_score": []}
    n_clusters_range = tuple(range(4, 100))
    for n_clusters in n_clusters_range:
        kmeans_model = KMeans(n_clusters=n_clusters, random_state=RANDOM_STATE).fit(
            embeddings
        )
        silhouette_avg = silhouette_score(embeddings, kmeans_model.labels_)
        results["n_clusters"].append(n_clusters)
        results["silhouette_score"].append(float(silhouette_avg))
        print(
            f"n_clusters = {n_clusters} average silhouette_score is : {silhouette_avg}"
        )

    to_yaml(results, "../data/kmeans-scan.yaml")


if __name__ == "__main__":
    main()
