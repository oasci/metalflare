#!/usr/bin/env python3
import json
import os
import string

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Input files
data_paths = {
    "Reduced": "../data/reduced-vamp.npy",
    "Oxidized": "../data/oxidized-vamp.npy",
    "Cu": "../data/cu-vamp.npy",
}

# Settings
n_clusters = 11
output_json = "../data/vamp-kmeans-global.json"
output_fig = "../data/vamp-kmeans-global.png"
output_label_template = "../data/{}-kmeans-labels.npy"

# === LOAD AND CONCATENATE ALL DATA === #
data_concat = []
individual_data = {}
for name, path in data_paths.items():
    data = np.load(path, mmap_mode="r")[:, :2]
    individual_data[name] = data
    data_concat.append(data)

data_all = np.vstack(data_concat)

data_clustering = data_all

# === FIT KMEANS === #
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
kmeans.fit(data_clustering)


def sort_centers_spatially(centers):
    centers = np.array(centers)
    n = len(centers)
    visited = []
    remaining = list(range(n))

    # Start at the point with the lowest x + y
    current = min(remaining, key=lambda i: centers[i][0] + centers[i][1])
    visited.append(current)
    remaining.remove(current)

    while remaining:
        dists = cdist([centers[current]], centers[remaining])
        next_idx = remaining[np.argmin(dists)]
        visited.append(next_idx)
        remaining.remove(next_idx)
        current = next_idx

    return visited


# === SORT CLUSTER CENTERS SPATIALLY === #
sorted_indices = sort_centers_spatially(kmeans.cluster_centers_)
letter_labels = list(string.ascii_uppercase)[:n_clusters]

cluster_centers = {
    letter_labels[i]: kmeans.cluster_centers_[idx].tolist()
    for i, idx in enumerate(sorted_indices)
}
label_map = {idx: letter_labels[i] for i, idx in enumerate(sorted_indices)}

# === SAVE TO JSON === #
with open(output_json, "w") as f:
    json.dump({"cluster_centers": cluster_centers}, f, indent=2)
print(f"Saved cluster center info to {output_json}")

# === BUILD CLUSTER REMAP: original idx -> sorted spatial idx === #
remap = {
    original_idx: sorted_idx for sorted_idx, original_idx in enumerate(sorted_indices)
}

# === SAVE REMAPPED LABELS === #
for name, data in individual_data.items():
    raw_labels = kmeans.predict(data)  # Raw KMeans labels (0 to n_clusters-1)

    # Remap using spatial sort: e.g., cluster originally idx=3, now becomes 0 if it's the "A" cluster
    remapped_labels = np.array([remap[label] for label in raw_labels], dtype=int)

    # Save remapped numeric labels
    output_path = output_label_template.format(name.lower().replace(" ", ""))
    np.save(output_path, remapped_labels)
    print(f"Saved remapped cluster labels for {name} to {output_path}")

# === PREDICT ALL LABELS === #
labels_all = kmeans.predict(data_all)
label_letters_all = [label_map[l] for l in labels_all]

# === QUICK PLOT === #
fig, ax = plt.subplots(figsize=(6, 5))
scatter = ax.scatter(
    data_all[:, 0], data_all[:, 1], c=np.array(labels_all), s=5, alpha=0.5, cmap="tab10"
)

# Annotate cluster centers with letters
for letter, coords in cluster_centers.items():
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

ax.set_xlabel("VAMP 1")
ax.set_ylabel("VAMP 2")
ax.grid(True, linestyle="--", alpha=0.3)
fig.tight_layout()
fig.savefig(output_fig, dpi=300)
print(f"Saved preview figure to {output_fig}")
