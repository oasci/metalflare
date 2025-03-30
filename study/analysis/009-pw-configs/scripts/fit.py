#!/usr/bin/env python3
import os

import numpy as np
import pandas as pd
from deeptime.decomposition import VAMP

os.chdir(os.path.dirname(os.path.realpath(__file__)))

base_dir = "../../../"

path_data = "../data/feats-proton-wire.parquet"

states = {
    0: "reduced",
    1: "oxidized",
    2: "cu",
}

df = pd.read_parquet(path_data)

X = df.to_numpy()[:, 1:]

vamp = VAMP(lagtime=1, dim=2)

model = vamp.fit_from_timeseries(X).fetch_model()

projected = model.transform(X)

for state, label in states.items():
    idx_state = df.index[df["state"] == state].to_numpy()
    data_state = projected[idx_state, :]
    np.save(f"../data/{label}-vamp.npy", data_state)
