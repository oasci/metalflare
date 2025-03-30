#!/usr/bin/env python3
import json
import os

import numpy as np
from scipy.spatial.distance import cdist

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# === CONFIGURATION === #
condition_labels = {"Reduced": "reduced", "Oxidized": "oxidized", "Cu(I)": "cu"}

# Parameters
base_dir = "../../../"
cluster_json_path = "../data/vamp-kmeans-global.json"
data_paths = {
    "Reduced": "../data/reduced-vamp.npy",
    "Oxidized": "../data/oxidized-vamp.npy",
    "Cu(I)": "../data/cu-vamp.npy",
}
output_json_path = "../data/cluster-center-closest-frames.json"

density_radius = 0.2  # cutoff for local density

# === LOAD CLUSTER CENTERS === #
with open(cluster_json_path, "r") as f:
    cluster_centers = json.load(f)["cluster_centers"]

cluster_coords = {label: np.array(coords) for label, coords in cluster_centers.items()}

# === LOAD ALL DATA === #
all_data = {}
for label, path in data_paths.items():
    all_data[label] = np.load(path, mmap_mode="r")

# === ANALYZE EACH CLUSTER CENTER === #
results = {}

for cluster_label, center in cluster_coords.items():
    condition_densities = {}
    best_score = -np.inf
    best_condition = None
    best_index = None

    for condition, data in all_data.items():
        distances = cdist([center], data)[0]
        local_count = np.sum(distances < density_radius)
        total_frames = len(data)
        percent_density = local_count / total_frames
        condition_key = condition_labels[condition]
        condition_densities[condition_key] = round(percent_density, 5)

        closest_idx = np.argmin(distances)
        closest_distance = distances[closest_idx]

        score = percent_density - closest_distance * 10.0

        if score > best_score:
            best_score = score
            best_condition = condition
            best_index = int(closest_idx)

    results[cluster_label] = {
        "condition": condition_labels[best_condition],  # type: ignore
        "frame": best_index,
        "densities": condition_densities,
    }

# === SAVE RESULTS === #
with open(output_json_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"Saved cluster-frame mapping with full density info to: {output_json_path}")
