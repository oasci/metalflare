#!/usr/bin/env python3
from typing import Any

import json
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from deeptime.covariance import Covariance
from deeptime.decomposition import VAMP
from sklearn.model_selection import TimeSeriesSplit

from metalflare.analysis.figures import use_mpl_rc_params

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Specify the paths to the trajectory and topology files
base_dir = "../../../"

# Paths
path_data = "../data/feats-proton-wire.parquet"
output_file = "../data/vamp_e_cv_results.json"

# Parameters and settings
lags_to_test = list(range(2, 151))
n_splits = 5  # number of folds for TimeSeriesSplit
dim_latent = 2  # Number of latent dimensions
t_delta = 10  # ps per time step

# Load data
df = pd.read_parquet(path_data)
# Each element in Xs is a numpy array for a unique state-run combination (features only).
# Trims the state and run labels
Xs = [group.to_numpy()[:, 2:] for _, group in df.groupby(["state", "run"])]

# Fitting loop
results: dict[int | str, Any] = {}
# Evaluate each lag value via cross-validation
for lag in lags_to_test:
    scores_fold = []

    # For each fold index, collect train and test segments across all runs.
    for fold in range(n_splits):
        Xs_train = []
        Xs_test = []

        # Process each run individually
        for run in Xs:
            tscv = TimeSeriesSplit(n_splits=n_splits)
            splits = list(tscv.split(run))
            # Ensure the run has enough samples for the current fold
            if fold < len(splits):
                train_idx, test_idx = splits[fold]
                Xs_train.append(run[train_idx])
                Xs_test.append(run[test_idx])

        # Only if we collected training and testing segments from at least one run
        if Xs_train and Xs_test:
            # Compute variances for train and test VAMP
            cov_train = Covariance(
                lagtime=lag, compute_c00=True, compute_c0t=True, compute_ctt=True
            ).fit_fetch(Xs_train)
            cov_test = Covariance(
                lagtime=lag, compute_c00=True, compute_c0t=True, compute_ctt=True
            ).fit_fetch(Xs_test)

            # Fit one VAMP model on the combined training data
            vamp_train = (
                VAMP(lagtime=lag, dim=dim_latent)
                .fit_from_covariances(cov_train)
                .fetch_model()
            )
            # Similarly, fit one model on the test data for scoring
            vamp_test = (
                VAMP(lagtime=lag, dim=dim_latent)
                .fit_from_covariances(cov_test)
                .fetch_model()
            )

            # Compute VAMP-E score for this fold
            score = vamp_train.score(r="VAMPE", test_model=vamp_test)
            scores_fold.append(score)

    # Aggregate scores for this lag across all folds
    results[lag] = {
        "vamp_e_scores": scores_fold,
        "mean_score": float(np.mean(scores_fold)) if scores_fold else None,
        "std_score": float(np.std(scores_fold)) if scores_fold else None,
    }

# Select the optimal lag based on the highest mean VAMP-E score
# (ignoring any lag with no valid scores)
valid_lags = {lag: res for lag, res in results.items() if res["mean_score"] is not None}
best_lag = int(max(valid_lags, key=lambda l: valid_lags[l]["mean_score"]))
results["best_lag"] = best_lag
print(f"Optimal lag: {best_lag} ({(t_delta*best_lag)} ps)")

# Fit the final model on all the data using the best lag
cov_all = Covariance(
    lagtime=best_lag, compute_c00=True, compute_c0t=True, compute_ctt=True
).fit_fetch(Xs)
best_model = (
    VAMP(lagtime=best_lag, dim=dim_latent).fit_from_covariances(cov_all).fetch_model()
)
projected = best_model.transform(df.to_numpy()[:, 2:])

# Save individual VAMP latent spaces based on state
states = {
    0: "reduced",
    1: "oxidized",
    2: "cu",
}
for state, label in states.items():
    idx_state = df.index[df["state"] == state].to_numpy()
    data_state = projected[idx_state, :]
    np.save(f"../data/{label}-vamp.npy", data_state)

# Save results
with open(output_file, "w") as f:
    json.dump(results, f, indent=2)
print(f"Cross-validated VAMP-E results saved to: {output_file}")

# Save the best model parameters
best_model_params = best_model.get_params()
params_file = "../data/best_model_params.pkl"
with open(params_file, "wb") as f:
    pickle.dump(best_model_params, f)
print(f"Best model parameters saved to: {params_file}")


# Update plot params
rc_json_path = os.path.join(base_dir, "misc/003-figure-style/matplotlib-rc-params.json")
font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
use_mpl_rc_params(rc_json_path, font_dirs)  # type: ignore

# Plotting the mean VAMP-E scores vs. lag values
# Extract lag values and their corresponding mean scores.
lags = sorted([lag for lag in results if isinstance(lag, int)])
mean_scores = [results[lag]["mean_score"] for lag in lags]

fig, ax = plt.subplots(1, 1, figsize=(3.5, 3))
ax.plot([lag * t_delta for lag in lags], mean_scores, marker="", linestyle="-")
ax.set_xlabel("Lag Time (ps)")
ax.set_ylabel("Mean VAMP-E Score")
ax.set_xlim(float(min(lags_to_test) * t_delta), float(max(lags_to_test) * t_delta))
plt.tight_layout()

# Save the plot in the data directory.
plot_path = "../data/vamp_e_cv_plot.png"
plt.savefig(plot_path)
plt.close()
print(f"Plot saved to: {plot_path}")
