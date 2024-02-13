#!/usr/bin/env python3

import os
import itertools

import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.realpath(__file__)))

features_path = "../../989-rogfp-cu-md-clustering/data/sim-features.parquet"

target_features = {
    "thr201_o-cym202_sg-dist": [3.15, 3.56, None],
    "thr201_hg1_og1_cb_cg2-dihedral": [-124.86, 44.50, 116.58, None],
    "ser203_og_cb_ca_n-dihedral": [-169.55, -73.33, -5.95, 69.73, None],
    "ser203_hg_og_cb_ca-dihedral": [-132.79, -82.7, -66.13, 58.20, 135.32, None],
    "cro65_og1_cb1_ca1_c1-dihedral": [-174.95, -68.29, 48.11, None],
    "cym202_c-ser203_n_ca_cb-dihedral": [-178.56, -79.82, -71.89, 126.31, None],
    "ser203_h-asn144_o-dist": [2.07, 5.24, 6.16, None],
}

def feat_angles(df):
    for col in df.columns:
        if "dihedral" in col:
            df[col] = np.degrees(df[col])
    return df


def find_representatives(features, target_features):
    # Filter out None values and get the mask of valid target features
    valid_mask = np.array([tf is not None for tf in target_features])
    valid_features = features[:, valid_mask]

    valid_target_features = np.array([tf for tf in target_features if tf is not None])

    # Calculate the Euclidean distances from each row in the dataset to the target features
    distances = np.linalg.norm(valid_features - valid_target_features, axis=1)

    # Find the index of the row with the minimum distance
    ranked_idxs = np.argsort(distances)
    distances = distances[ranked_idxs]

    return ranked_idxs, distances



def main():
    df = pd.read_parquet(features_path)
    df = feat_angles(df)
    df_rogfp = df[df["system_label"] == "rogfp"]
    df_rogfp_cu = df[df["system_label"] == "rogfp_cu"]
    feat_rogfp = df_rogfp[target_features.keys()].to_numpy()
    feat_rogfp_cu = df_rogfp_cu[target_features.keys()].to_numpy()

    combinations = itertools.product(*target_features.values())

    structure_info = {feature: [] for feature in target_features.keys()}
    structure_info["system_label"] = []
    structure_info["distance"] = []
    structure_info["traj_idx"] = []
    for combination in combinations:
        # Skip if all target features are None
        if len(set(combination)) == 1:
            continue
        target_feature = np.array(combination)
        rogfp_rep_idxs, rogfp_rep_dist = find_representatives(feat_rogfp, target_feature)
        for key, value in zip(target_features.keys(), combination):
            structure_info[key].append(value)
        structure_info["system_label"].append("rogfp")
        structure_info["traj_idx"].append(rogfp_rep_idxs[0])
        structure_info["distance"].append(rogfp_rep_dist[0])

        rogfp_cu_rep_idxs, rogfp_cu_rep_dist = find_representatives(feat_rogfp_cu, target_feature)
        for key, value in zip(target_features.keys(), combination):
            structure_info[key].append(value)
        structure_info["system_label"].append("rogfp_cu")
        structure_info["traj_idx"].append(rogfp_cu_rep_idxs[0])
        structure_info["distance"].append(rogfp_cu_rep_dist[0])

    df_structure_info = pd.DataFrame(structure_info)
    df_structure_info.to_csv("../data/representative_structures.csv", index=False)


if __name__ == "__main__":
    main()
