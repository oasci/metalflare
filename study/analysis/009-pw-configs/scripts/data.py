#!/usr/bin/env python3
import os

import pandas as pd

from metalflare.utils import load_features

os.chdir(os.path.dirname(os.path.realpath(__file__)))

data_y_str = "cro65_oh-his146_hd1-dist"
data_y_label = "Cro66 OH - His148 HD1 Distance [Å]"

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}
names_data = [
    # Sensing Cys backbone
    # "asn144_c-cys145_n_ca_c-dihedral",
    # "cys145_n_ca_c-his146_n-dihedral",
    # "thr201_c-cys202_n_ca_c-dihedral",
    # "cys202_n_ca_c-ser203_n-dihedral",
    # Cro dihedral
    "cro65_og1_cb1_ca1_c1-dihedral",
    # Proton wire backbone
    "ser200_c-thr201_n_ca_c-dihedral",
    "thr201_n_ca_c-cys202_n-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
    "leu219_c-glu220_n_ca_c-dihedral",
    "glu220_n_ca_c-phe221_n-dihedral",
    # Proton wire distances
    "ser203_og-glu220_he2-dist",
    # To CRO distances
    # "cro65_oh-thr201_hg1-dist",
    # "cro65_oh-tyr143_hh-dist",
]

base_dir = "../../../"


dfs = []
idx = 0
for state_key, state_path in names_state.items():
    print(f"Working on {state_key}...")
    paths_data = [
        os.path.join(base_dir, f"analysis/{state_path}/data/struct-desc/{name}.npy")
        for name in names_data
    ]

    df = load_features(paths_data, normalize_distances=True, transform_dihedrals=True)
    df.insert(0, "state", idx)
    dfs.append(df)

    idx += 1

df_comb = pd.concat(dfs, ignore_index=True)

df_comb.to_parquet("../data/feats-proton-wire.parquet")
