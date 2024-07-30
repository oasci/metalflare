#!/usr/bin/env python3
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler, MinMaxScaler

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
    "asn142_c-tyr143_n_ca_c-dihedral",
    "tyr143_n_ca_c-asn144_n-dihedral",
    "asn144_c-cys145_n_ca_c-dihedral",
    "cys145_n_ca_c-his146_n-dihedral",
    "asn144_n_ca_c-cys145_n-dihedral",
    "tyr143_c-asn144_n_ca_c-dihedral",
    "cys145_c-his146_n_ca_c-dihedral",
    "his146_n_ca_c-asn147_n-dihedral",
    "cys202_c-ser203_n_ca_c-dihedral",
    "ser203_n_ca_c-ala204_n-dihedral",
    "cys202_n_ca_c-ser203_n-dihedral",
    "thr201_c-cys202_n_ca_c-dihedral",
    "ser200_c-thr201_n_ca_c-dihedral",
    "thr201_n_ca_c-cys202_n-dihedral",
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
    scaler = StandardScaler()
    normalized_df = pd.DataFrame(
        scaler.fit_transform(combined_df), columns=combined_df.columns
    )

    return normalized_df


def train_and_evaluate_model(model, X_train, y_train):
    model.fit(X_train, y_train)
    return model


def get_normalized_feature_importance(model, X, feature_names):
    if hasattr(model, 'feature_importances_'):
        importance = model.feature_importances_
    elif hasattr(model, 'coef_'):
        importance = np.abs(model.coef_)
    else:
        raise ValueError("Model does not have feature_importances_ or coef_ attribute")

    # Reshape importance to 2D array for MinMaxScaler
    importance_2d = importance.reshape(-1, 1)

    # Initialize and fit MinMaxScaler
    scaler = MinMaxScaler(feature_range=(0, 1))
    normalized_importance = scaler.fit_transform(importance_2d).flatten()

    return pd.DataFrame(
        {"Feature": feature_names, "Importance": normalized_importance}
    ).sort_values("Importance", ascending=False)


def generate_markdown_table(df):
    return df.to_markdown(index=False, floatfmt=".4f")


if __name__ == "__main__":
    base_dir = "../../../"

    for state_key, state_path in names_state.items():
        print(f"Working on {state_key}...")
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

        # Make dataset smaller
        random_seed = 5478420
        # X, _, y, _ = train_test_split(X, y, test_size=0.1)

        models = {
            "Random Forest": RandomForestRegressor(
                n_estimators=100, random_state=random_seed
            ),
            "Gradient Boosting": GradientBoostingRegressor(
                n_estimators=100, random_state=random_seed
            ),
            "Ridge": Ridge(alpha=1.0),
        }

        importance_dfs = {}
        for model_name, model in models.items():
            trained_model = train_and_evaluate_model(model, X, y)
            importance_df = get_normalized_feature_importance(
                trained_model, X, feature_names
            )
            importance_df.columns = ["Feature", model_name]
            importance_dfs[model_name] = importance_df

        # Combine all importance DataFrames
        combined_importance = (
            importance_dfs["Random Forest"]
            .merge(importance_dfs["Gradient Boosting"], on="Feature")
            .merge(importance_dfs["Ridge"], on="Feature")
        )

        # Calculate average importance and sort
        combined_importance["Mean Importance"] = combined_importance[
            ["Random Forest", "Gradient Boosting", "Ridge"]
        ].mean(axis=1)
        combined_importance = combined_importance.sort_values(
            "Mean Importance", ascending=False
        )

        # Reorder columns
        column_order = [
            "Feature",
            "Mean Importance",
            "Random Forest",
            "Gradient Boosting",
            "Ridge",
        ]
        combined_importance = combined_importance[column_order]

        # Generate report
        report = f"{generate_markdown_table(combined_importance)}\n"

        # Write the report to a file
        with open(f"{state_key}-feature-report.md", "w") as f:
            f.write(report)
