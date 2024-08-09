#!/usr/bin/env python3
import json
import os

import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

from metalflare.analysis.figures import use_mpl_rc_params
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
    # Residue specific
    # Relevant backbones
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
    "glu220_n_ca_c-phe221_n-dihedral",
]


def optimize_model(X, y, model, param_grid, cv=3):
    pipeline = Pipeline([("scaler", StandardScaler()), ("model", model)])

    grid_search = GridSearchCV(
        pipeline,
        param_grid,
        cv=cv,
        n_jobs=-1,
        verbose=1,
        scoring="neg_mean_squared_error",
    )
    grid_search.fit(X, y)

    return grid_search.best_estimator_, grid_search.best_params_


if __name__ == "__main__":
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    random_seed = 3984729

    # Define models and their parameter grids
    models = {
        "XGBoost": (
            XGBRegressor(random_state=random_seed),
            {
                "model__n_estimators": [250, 400, 550, 700],
                "model__learning_rate": [0.05, 0.1, 0.2],
                "model__max_depth": [5, 7, 9],
                "model__reg_alpha": [0.0, 0.1, 0.2],
                "model__reg_lambda": [1.0, 0.9, 0.8],
            },
        ),
        "ElasticNet": (
            ElasticNet(random_state=random_seed),
            {
                "model__alpha": [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 2.0, 5.0],
                "model__l1_ratio": [0.2, 0.5, 0.8, 1.0],
                "model__max_iter": [100000],
            },
        ),
    }

    all_results = {}

    for state_key, state_path in names_state.items():
        print(f"\nProcessing state: {state_key}")

        # Load data for the current state
        paths_data = [
            os.path.join(base_dir, f"analysis/{state_path}/data/struct-desc/{name}.npy")
            for name in names_data
        ]
        X = load_features(paths_data, transform_dihedrals=True)
        y = np.load(
            os.path.join(
                base_dir, f"analysis/{state_path}/data/struct-desc/{data_y_str}.npy"
            )
        )

        # Split the data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=random_seed
        )

        state_results = {}

        # Optimize models, evaluate performance, and calculate feature importance
        for model_name, (model, param_grid) in models.items():
            print(f"Optimizing {model_name} for {state_key}...")
            best_model, best_params = optimize_model(
                X_train, y_train, model, param_grid
            )

            # Evaluate the model
            y_pred = best_model.predict(X_test)
            mse = mean_squared_error(y_test, y_pred)
            r2 = r2_score(y_test, y_pred)

            # Get feature importance
            if model_name == "XGBoost":
                importance = best_model.named_steps["model"].feature_importances_
            elif model_name == "ElasticNet":
                importance = np.abs(best_model.named_steps["model"].coef_)

            state_results[model_name] = {
                "best_params": best_params,
                "mse": mse,
                "r2": r2,
                "feature_importance": dict(zip(X.columns, importance.tolist())),
            }

            print(f"Best parameters for {model_name}: {best_params}")
            print(f"MSE: {mse:.4f}")
            print(f"R2 Score: {r2:.4f}")

        all_results[state_key] = state_results

    # Save results to JSON file
    output_file = "model_results_by_state.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=4)

    print(f"\nResults saved to {output_file}")
