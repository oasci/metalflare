#!/usr/bin/env python3
import os

import pandas as pd
import numpy as np

from metalflare.utils import load_features

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Each simulation states has three concatenated trajectories that are the same
# length, but with different initial conditions.
# We need to mark these trajectories boundaries to avoid any nonphysical
# transitions.
n_independent_runs = 3

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}
names_data = [
    ###   DIHEDRALS   ###
    # Cro
    "cro65_og1_cb1_ca1_c1-dihedral",
    # Proton wire backbone
    "ser200_c-thr201_n_ca_c-dihedral",
    "thr201_n_ca_c-cys202_n-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
    "leu219_c-glu220_n_ca_c-dihedral",
    "glu220_n_ca_c-phe221_n-dihedral",
    # glu222
    "glu220_ca_cb_cg_cd-dihedral",
    "glu220_cb_cg_cd_oe2-dihedral",
    # ser205
    "ser203_n_ca_cb_og-dihedral",
    # thr203
    "thr201_n_ca_cb_og1-dihedral",
    ###   DISTANCES   ###
    # To CRO
    "cro65_oh-thr201_hg1-dist",
    "cro65_og1-glu220_he2-dist",
    # Proton wire
    "ser203_og-glu220_he2-dist",
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
    n_steps_per_traj = int(df.shape[0] / n_independent_runs)
    df.insert(1, "run", np.repeat(np.arange(n_independent_runs), n_steps_per_traj))
    dfs.append(df)

    idx += 1

df_comb = pd.concat(dfs, ignore_index=True)

df_comb.to_parquet("../data/feats-proton-wire.parquet")
