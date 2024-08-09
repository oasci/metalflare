#!/usr/bin/env python3
import os

import numpy as np

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pls import (
    compare_states,
    compute_loading_angles_and_magnitudes,
    create_angle_magnitude_file,
    perform_pls_regression,
    plot_pls_results,
)
from metalflare.utils import load_features

os.chdir(os.path.dirname(os.path.realpath(__file__)))

data_y_str = "cro65_oh-tyr143_hh-dist"
data_y_label = "Cro66 OH - Tyr145 HH Distance [Å]"

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}
names_data = [
    # Residue specific
    "tyr143_ca_cb_cg_cd1-dihedral",
    "tyr143_ce1_cz_oh_hh-dihedral",
    "tyr143_hh-thr62_og1-dist",
    # Backbone distances
    "his146_h-thr201_o-dist",
    "ser203_h-asn144_o-dist",
    "cys202_o-phe221_h-dist",
    # Backbone angles
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
    "cys202_n_ca_c-ser203_n-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
    "leu219_c-glu220_n_ca_c-dihedral",
    "glu220_n_ca_c-phe221_n-dihedral"
    # Backbone distances
    # "cys202_o-phe221_h-dist",
    # "his146_h-thr201_o-dist",
    # "ser203_h-asn144_o-dist",
    # Non-backbone distances
    # "ser203_og-glu220_he2-dist",
    # "his146_hd1-asn144_o-dist",
    # "cro65_og1-glu220_he2-dist",
]

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

        X = load_features(paths_data, normalize_distances=True, transform_dihedrals=True)
        y = np.load(
            os.path.join(
                base_dir, f"analysis/{state_path}/data/struct-desc/{data_y_str}.npy"
            )
        )
        feature_names = X.columns

        pls, X_scaled, y_scaled, scaler_X, scaler_y = perform_pls_regression(X, y)
        gradient_vector = plot_pls_results(
            pls,
            X_scaled,
            y_scaled,
            feature_names,
            state_key,
            y_label=data_y_label,
            arrow_scaling=4.0,
        )

        if gradient_vector is not None:
            angle_magnitude_data = compute_loading_angles_and_magnitudes(
                pls, gradient_vector, feature_names
            )
            create_angle_magnitude_file(angle_magnitude_data, state_key, feature_names)

        results[state_key] = (pls, X_scaled, y_scaled, scaler_X, scaler_y)

    compare_states(results, feature_names)
