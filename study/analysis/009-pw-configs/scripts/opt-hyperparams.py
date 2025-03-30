#!/usr/bin/env python3
import json
import os

import numpy as np
import pandas as pd
from deeptime.decomposition import VAMP
from sklearn.model_selection import KFold

# Setup
os.chdir(os.path.dirname(os.path.realpath(__file__)))
path_data = "../data/feats-proton-wire.parquet"

# Parameters
lags_to_test = list(range(1, 21))
lags_to_test.extend(list(range(25, 51, 5)))
n_splits = 5
dim = 2
output_file = "../data/vamp_e_cv_results.json"

# Load data
df = pd.read_parquet(path_data)
X = df.to_numpy()[:, 1:]

results = {}
kf = KFold(n_splits=n_splits, shuffle=False)

for lag in lags_to_test:
    scores = []

    for train_idx, test_idx in kf.split(X):
        X_train = X[train_idx]
        X_test = X[test_idx]

        # Fit on training data
        vamp_train = (
            VAMP(lagtime=lag, dim=dim).fit_from_timeseries(X_train).fetch_model()
        )

        # Fit on test data
        vamp_test = VAMP(lagtime=lag, dim=dim).fit_from_timeseries(X_test).fetch_model()

        # Compute cross-validated VAMP-E score
        score = vamp_train.score(r="VAMPE", test_model=vamp_test)
        scores.append(score)

    results[lag] = {
        "vamp_e_scores": scores,
        "mean_score": float(np.mean(scores)),
        "std_score": float(np.std(scores)),
    }

# Select optimal lag
best_lag = max(results, key=lambda l: results[l]["mean_score"])
results["best_lag"] = best_lag

# Save results
with open(output_file, "w") as f:
    json.dump(results, f, indent=2)

print(f"Cross-validated VAMP-E results saved to: {output_file}")
print(f"Optimal lag (based on VAMP-E): {best_lag}")
