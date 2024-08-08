#!/usr/bin/env python3
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.utils import format_feature_name

os.chdir(os.path.dirname(os.path.realpath(__file__)))

results_paths = {
    "all": "./model_results.json",
    "individual": "./model_results_by_state.json",
}


def plot_feature_importance(importance, names, model_type, system):
    feature_importance = np.array(importance)
    feature_names = np.array(names)

    data = {"feature_names": feature_names, "feature_importance": feature_importance}
    fi_df = pd.DataFrame(data)

    fi_df.sort_values(by=["feature_importance"], ascending=False, inplace=True)

    plt.figure(figsize=(6, 5))
    plt.bar(range(fi_df.shape[0]), fi_df["feature_importance"])
    plt.xticks(range(fi_df.shape[0]), fi_df["feature_names"], rotation="vertical")
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.tight_layout()
    plt.savefig(f"{system}_{model_type.lower()}_feature_importance.png")
    plt.close()


# In the main part of the script, modify the call to compare_states:
if __name__ == "__main__":
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    for kind, json_path in results_paths.items():
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if kind == "all":
            data = {"all": data}

        for system, model_info in data.items():
            for model_name, results in model_info.items():
                feature_names = list(results["feature_importance"].keys())
                feature_names = [
                    format_feature_name(i).replace(" dihedral", "") for i in feature_names
                ]

                plot_feature_importance(
                    list(results["feature_importance"].values()),
                    feature_names,
                    model_name,
                    system,
                )
