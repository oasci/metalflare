#!/usr/bin/env python3
import json
import os
from itertools import combinations

import numpy as np

from metalflare.utils import format_feature_name, load_features


def get_correlation(x, y):
    return np.corrcoef(x, y)[0, 1]


os.chdir(os.path.dirname(os.path.realpath(__file__)))

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}

names_data = [
    "tyr143_ca_cb_cg_cd1-dihedral",
    "tyr143_ce1_cz_oh_hh-dihedral",
    "asn142_c-tyr143_n_ca_c-dihedral",
    "tyr143_n_ca_c-asn144_n-dihedral",
    "tyr143_c-asn144_n_ca_c-dihedral",
    "asn144_n_ca_c-cys145_n-dihedral",
    "asn144_c-cys145_n_ca_c-dihedral",
    "cys145_n_ca_c-his146_n-dihedral",
    "cys145_c-his146_n_ca_c-dihedral",
    "his146_n_ca_c-asn147_n-dihedral",
    "ser200_c-thr201_n_ca_c-dihedral",
    "thr201_n_ca_c-cys202_n-dihedral",
    "thr201_c-cys202_n_ca_c-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
]


def process_trajectory(X):
    traj_results = {}
    for (i, name1), (j, name2) in combinations(enumerate(names_data), 2):
        label = f"{format_feature_name(name1)} to {format_feature_name(name2)}"
        corr = get_correlation(X.iloc[:, i].to_numpy(), X.iloc[:, j].to_numpy())
        traj_results[label] = corr
    return traj_results


def print_statistics(data):
    all_correlations = []
    for state in data.values():
        all_correlations.extend(list(state.values()))

    all_correlations = np.abs(all_correlations)

    mean = np.mean(all_correlations)
    median = np.median(all_correlations)
    std_dev = np.std(all_correlations)

    cutoff_1std = mean + std_dev
    cutoff_2std = mean + 2 * std_dev
    cutoff_3std = mean + 3 * std_dev

    print("\nStatistics:")
    print(f"Mean absolute correlation: {mean:.4f}")
    print(f"Median absolute correlation: {median:.4f}")
    print(f"Standard Deviation: {std_dev:.4f}")
    print(f"\nPotential cutoffs:")
    print(f"  1 std dev above mean: {cutoff_1std:.4f}")
    print(f"  2 std dev above mean: {cutoff_2std:.4f}")
    print(f"  3 std dev above mean: {cutoff_3std:.4f}")


if __name__ == "__main__":
    base_dir = "../../../"
    all_results = {}
    for state_key, state_path in names_state.items():
        print(f"\nProcessing state: {state_key}")
        state_results = {}
        paths_data = [
            os.path.join(base_dir, f"analysis/{state_path}/data/struct-desc/{name}.npy")
            for name in names_data
        ]
        X_full = load_features(paths_data, transform_dihedrals=True)

        # Split the concatenated data into three trajectories
        traj_length = len(X_full) // 3
        X_trajectories = [
            X_full.iloc[i * traj_length : (i + 1) * traj_length] for i in range(3)
        ]

        for traj_num, X in enumerate(X_trajectories, 1):
            print(f"Processing trajectory {traj_num}")
            traj_results = process_trajectory(X)

            # Accumulate results
            for pair, corr in traj_results.items():
                if pair not in state_results:
                    state_results[pair] = []
                state_results[pair].append(corr)

        # Average correlations across trajectories
        for pair, corrs in state_results.items():
            state_results[pair] = np.mean(corrs, axis=0)

        all_results[state_key] = state_results

    # Prepare data for JSON
    data = {}
    for state, correlations in all_results.items():
        data[state] = {}
        # Sort correlations by absolute value in descending order
        sorted_correlations = sorted(
            correlations.items(), key=lambda x: abs(x[1]), reverse=True
        )
        for pair, value in sorted_correlations:
            data[state][pair] = value

    # Print statistics
    print_statistics(data)

    # Save results to JSON file
    with open("dihedral-correlation.json", "w+") as f:
        json.dump(data, f, indent=4)

print("\nResults saved to dihedral-correlation.json")
