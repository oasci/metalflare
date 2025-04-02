#!/usr/bin/env python3
import json
import os
import string

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from sklearn import preprocessing

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Specify the paths to the trajectory and topology files
base_dir = "../../../"

# Update plot params
rc_json_path = os.path.join(base_dir, "misc/003-figure-style/matplotlib-rc-params.json")
font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
use_mpl_rc_params(rc_json_path, font_dirs)  # type: ignore

# Input files
data_paths = {
    "Reduced": "../data/reduced-vamp.npy",
    "Oxidized": "../data/oxidized-vamp.npy",
    "Cu": "../data/cu-vamp.npy",
}

coords_clustering = np.array(
    [
        [-1.07, 3.04],  # A
        [-1.045, 1.9],  # B
        [-1.025, 1.27],  # C
        [-1.04, 0.9],  # D
        [-0.99, 0.1],  # E
        [-1.00, -0.74],  # F
        [-0.90, 1.1],  # G
    ]
)

# Settings
output_json = "../data/vamp-clustering-global.json"
output_fig = "../data/vamp-clustering-global.png"
output_label_template = "../data/{}-clustering-labels.npy"


letter_labels = list(string.ascii_uppercase)[: coords_clustering.shape[0]]

cluster_centers = {
    letter_labels[i]: coords_clustering[i].tolist()
    for i, idx in enumerate(letter_labels)
}
label_map = {i: letter_labels[i] for i, idx in enumerate(letter_labels)}

with open(output_json, "w") as f:
    json.dump({"cluster_centers": cluster_centers}, f, indent=2)
print(f"Saved cluster center info to {output_json}")


# === LOAD AND CONCATENATE ALL DATA === #
data_concat = []
individual_data = {}
for name, path in data_paths.items():
    data = np.load(path, mmap_mode="r")[:, :2]
    individual_data[name] = data
    data_concat.append(data)

data_all = np.vstack(data_concat)

scaler = preprocessing.MinMaxScaler().fit(data_all)
data_clustering = scaler.transform(data_all)
coords_clustering = scaler.transform(coords_clustering)

cluster_centers_normed = {
    letter_labels[i]: coords_clustering[i].tolist()
    for i, idx in enumerate(letter_labels)
}


def get_labels(data, coords_centroids):
    dists = cdist(data, coords_centroids)  # shape (n_samples, n_clusters)
    labels = np.argmin(dists, axis=1)
    return labels


for name, data in individual_data.items():
    data = scaler.transform(data)
    labels = get_labels(
        data, coords_clustering
    )  # Raw KMeans labels (0 to n_clusters-1)

    # Save remapped numeric labels
    output_path = output_label_template.format(name.lower().replace(" ", ""))
    np.save(output_path, labels)
    print(f"Saved remapped cluster labels for {name} to {output_path}")


labels_all = get_labels(data_clustering, coords_clustering)
label_letters_all = [label_map[l] for l in labels_all]


fig, ax = plt.subplots(figsize=(6, 5))
scatter = ax.scatter(
    data_clustering[:, 0],
    data_clustering[:, 1],
    c=np.array(labels_all),
    s=5,
    alpha=0.5,
    cmap="tab10",
)

# Annotate cluster centers with letters
for letter, coords in cluster_centers_normed.items():
    ax.text(
        coords[0],
        coords[1],
        letter,
        fontsize=12,
        weight="bold",
        ha="center",
        va="center",
        color="#ffffff",
    )

ax.set_xlabel("Normalized VAMP 1")
ax.set_ylabel("Normalized VAMP 2")
ax.grid(True, linestyle="--", alpha=0.3)
fig.tight_layout()
fig.savefig(output_fig, dpi=300)
print(f"Saved preview figure to {output_fig}")
