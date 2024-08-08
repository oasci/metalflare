#!/usr/bin/env python3
import json
import os
import matplotlib.pyplot as plt
import numpy as np
from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.utils import format_feature_name

os.chdir(os.path.dirname(os.path.realpath(__file__)))

results_paths = {
    "all": "./model_results.json",
    "individual": "./model_results_by_state.json",
}

def plot_combined_feature_importance(xgb_importance, en_importance, names, system):
    feature_names = np.array(names)
    xgb_importance = np.array(xgb_importance)
    en_importance = np.array(en_importance)

    # Sort features based on XGBoost importance
    sort_idx = np.argsort(xgb_importance)[::-1]
    feature_names = feature_names[sort_idx]
    xgb_importance = xgb_importance[sort_idx]
    en_importance = en_importance[sort_idx]

    x = np.arange(len(feature_names))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x - width/2, xgb_importance, width, label='XGBoost', color="#264653")
    ax.bar(x + width/2, en_importance, width, label='ElasticNet', color="#2a9d8f")

    ax.set_ylabel('Importance')
    ax.set_xlabel('Feature')
    ax.set_xticks(x)
    ax.set_xticklabels(feature_names, rotation='vertical')
    ax.set_xticklabels(feature_names, rotation=45, ha='right')
    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(f"{system}_feature_importance.png", dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    base_dir = "../../../"
    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    # Load the JSON data
    for kind, json_path in results_paths.items():
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if kind == "all":
            data = {"all": data}

        for system, model_info in data.items():
            xgb_importance = list(model_info["XGBoost"]["feature_importance"].values())
            en_importance = list(model_info["ElasticNet"]["feature_importance"].values())

            feature_names = list(model_info["XGBoost"]["feature_importance"].keys())
            feature_names = [
                format_feature_name(i).replace(" dihedral", "") for i in feature_names
            ]

            plot_combined_feature_importance(
                xgb_importance,
                en_importance,
                feature_names,
                system,
            )

    print("Feature importance plots have been generated for all systems.")
