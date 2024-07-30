#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import string
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

data_y_str = "cro65_oh-tyr143_hh-dist"
data_y_label = "Cro66 OH - Tyr145 HH Distance [Å]"

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}
names_data = [
    "tyr143_ca_cb_cg_cd1-dihedral",
    "tyr143_ce1_cz_oh_hh-dihedral",
    # "asn142_c-tyr143_n_ca_c-dihedral",
    "tyr143_n_ca_c-asn144_n-dihedral",
    "cys145_n_ca_c-his146_n-dihedral",
    "asn144_n_ca_c-cys145_n-dihedral",
    "tyr143_c-asn144_n_ca_c-dihedral",
    "cys145_c-his146_n_ca_c-dihedral",
    "his146_n_ca_c-asn147_n-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
    "ser200_c-thr201_n_ca_c-dihedral",

    "cro65_og1-glu220_he2-dist",
    # "cys202_o-phe221_h-dist",
    "his146_hd1-asn144_o-dist",
    "his146_h-thr201_o-dist",
    "ser203_h-asn144_o-dist",
    "ser203_og-glu220_he2-dist",
]

def load_and_preprocess_data(file_paths):
    df_list = []
    for path in file_paths:
        data = np.load(path)
        base_name = os.path.basename(path).split(".")[0]

        if "dihedral" in base_name:
            df = pd.DataFrame({f"{base_name}_cos": np.cos(data)})
        if "dist" in base_name:
            df = pd.DataFrame({base_name: data})

        df_list.append(df)

    combined_df = pd.concat(df_list, axis=1)
    return combined_df

def perform_pls_regression(X, y, n_components=2):
    scaler_X = StandardScaler()
    scaler_y = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)
    # y_scaled = scaler_y.fit_transform(y.reshape(-1, 1)).ravel()
    y_scaled = y

    pls = PLSRegression(n_components=n_components)
    pls.fit(X_scaled, y_scaled)

    return pls, X_scaled, y_scaled, scaler_X, scaler_y

def plot_pls_results(pls, X_scaled, y_scaled, feature_names, state_key):
    # Calculate PLS components
    pls_components = pls.transform(X_scaled)

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 10))
    scatter = ax.scatter(pls_components[:, 0], pls_components[:, 1], c=y_scaled, cmap='viridis', alpha=0.5)
    plt.colorbar(scatter, label=data_y_label)

    # Add labels and title
    ax.set_xlabel('PLS Component 1')
    ax.set_ylabel('PLS Component 2')

    # Create a dictionary to map features to letters
    feature_letters = dict(zip(feature_names, string.ascii_uppercase[:len(feature_names)]))

    offset_distance = 0.3

    for label, loading in zip(feature_names, pls.x_loadings_):
        letter = feature_letters[label]
        loading *= 4.0
        ax.arrow(
            0, 0, loading[0], loading[1],
            color='#ff686b', alpha=1.0, head_width=0.11, head_length=0.1,
            width=0.03
        )

        # Calculate unit vector in direction of arrow
        magnitude = np.sqrt(loading[0]**2 + loading[1]**2)
        unit_vector = loading / magnitude if magnitude != 0 else np.array([0, 0])

        # Calculate label position
        label_x = loading[0] + unit_vector[0] * offset_distance
        label_y = loading[1] + unit_vector[1] * offset_distance

        ax.annotate(letter, (label_x, label_y), xytext=(0, 0), textcoords='offset points',
            ha='center', va='center', alpha=1.0, fontweight='bold',
            bbox=dict(boxstyle='circle,pad=0.3', fc='white', ec='none', alpha=0.75)
        )

    # Add R-squared value
    r_squared = pls.score(X_scaled, y_scaled)
    ax.text(0.05, 0.95, f'R² = {r_squared:.3f}', transform=ax.transAxes, verticalalignment='top')

    plt.tight_layout()
    plt.savefig(f'{state_key}_pls_regression.png', dpi=300)
    plt.close()

    # Print the key to console
    print(f"\nKey for {state_key} PLS Regression Plot:")
    for feature, letter in feature_letters.items():
        print(f"-   **{letter}**: {feature}")

def compare_states(results, feature_names):
    # Extract loadings and create a DataFrame
    loadings_df = {}
    for state in results:
        loadings = results[state][0].x_loadings_
        # Calculate magnitude using the first two components
        magnitudes = np.sqrt(loadings[:, 0]**2 + loadings[:, 1]**2)
        loadings_df[state] = magnitudes

    loadings_df = pd.DataFrame(loadings_df)
    loadings_df.index = feature_names

    print("\nLargest loading magnitudes (based on first two components):")
    for state in ['reduced', 'oxidized', 'cu']:
        print(f"\n{state.capitalize()} state:")
        changes = loadings_df[state].sort_values(ascending=False)
        for feature, change in changes.head(10).items():
            print(f"{feature}: {change:.4f}")

    # # Calculate the difference in loadings
    # diff_df = loadings_df.sub(loadings_df['reduced'], axis=0)

    # # Print the most significant changes
    # print("\nMost significant changes in loading magnitudes (relative to reduced state):")
    # for state in ['oxidized', 'cu']:
    #     print(f"\n{state.capitalize()} state:")
    #     changes = diff_df[state].abs().sort_values(ascending=False)
    #     for feature, change in changes.head(5).items():
    #         print(f"{feature}: {change:.4f}")

# In the main part of the script, modify the call to compare_states:
if __name__ == "__main__":
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    results = {}

    for state_key, state_path in names_state.items():
        print(f"\nWorking on {state_key}...")
        paths_data = [
            os.path.join(base_dir, f"analysis/{state_path}/data/struct-desc/{name}.npy")
            for name in names_data
        ]

        X = load_and_preprocess_data(paths_data)
        y = np.load(
            os.path.join(
                base_dir, f"analysis/{state_path}/data/struct-desc/{data_y_str}.npy"
            )
        )
        feature_names = X.columns

        pls, X_scaled, y_scaled, scaler_X, scaler_y = perform_pls_regression(X, y)
        plot_pls_results(pls, X_scaled, y_scaled, feature_names, state_key)

        results[state_key] = (pls, X_scaled, y_scaled, scaler_X, scaler_y)

    compare_states(results, feature_names)
