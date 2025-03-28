#!/usr/bin/env python3
import os

import pandas as pd

from deeptime.decomposition import VAMP

os.chdir(os.path.dirname(os.path.realpath(__file__)))

data_y_str = "cro65_oh-his146_hd1-dist"
data_y_label = "Cro66 OH - His148 HD1 Distance [Å]"

names_state = {
    "reduced": "005-rogfp-glh-md",
    "oxidized": "007-rogfp-oxd-glh-md",
    "cu": "006-rogfp-cu-glh-md",
}


base_dir = "../../../"

path_data = "../data/feats-proton-wire.parquet"


df = pd.read_parquet(path_data)

X = df.to_numpy()[:, 1:]


vamp = VAMP(dim=1)

covars = vamp.covariance_estimator(lagtime=10).fit(X).fetch_model()

model = vamp.fit(covars).fetch_model()
