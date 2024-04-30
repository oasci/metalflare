#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

os.chdir(os.path.dirname(os.path.realpath(__file__)))

rogfp2_file_path = "../data/rogfp2-desc.parquet"
rogfp2_oxd_file_path = "../data/rogfp2-oxd-desc.parquet"
rogfp2_cu_file_path = "../data/rogfp2-cu-desc.parquet"


def feat_angles(df):
    """
    For each column in the DataFrame with 'dihedral' in its name, this function updates the column values from
    the range [-pi, pi] to [0, 2*pi] and replaces it with two new columns for sin(x) and cos(x) of the updated values.

    Parameters:
    - df: pandas.DataFrame containing at least one column with 'dihedral' in its name.

    Returns:
    - A modified DataFrame with the original 'dihedral' columns updated to [0, 2*pi] range and two new columns for
      each 'dihedral' column representing sin(x) and cos(x) of these values.
    """
    for col in df.columns:
        if "dihedral" in col:
            # Update the range from [-pi, pi] to [0, 2*pi]
            updated_values = df[col] + np.pi

            # Calculate sin(x) and cos(x)
            df[f"{col}_sin"] = np.sin(updated_values)
            df[f"{col}_cos"] = np.cos(updated_values)

            # Optional: Remove the original 'dihedral' column if not needed
            # df.drop(col, axis=1, inplace=True)

    return df


def featurization(df_data):
    # Transform angles to sin(x) and cos(x)
    df_data = feat_angles(df_data)

    # Scale features
    scaler = MinMaxScaler()
    data = scaler.fit_transform(df_data.values)

    return data


def main():
    df_rogfp2 = pd.read_parquet(rogfp2_file_path)
    df_rogfp2 = feat_angles(df_rogfp2)
    df_rogfp2["system_label"] = "rogfp_red"  # Add system label after processing

    df_oxd_rogfp2 = pd.read_parquet(rogfp2_oxd_file_path)
    df_oxd_rogfp2 = feat_angles(df_oxd_rogfp2)
    df_oxd_rogfp2["system_label"] = "rogfp_oxd"  # Add system label after processing

    df_rogfp2_cu = pd.read_parquet(rogfp2_cu_file_path)
    df_rogfp2_cu = feat_angles(df_rogfp2_cu)
    df_rogfp2_cu["system_label"] = "rogfp_cu"  # Add system label after processing

    # Combine the dataframes
    df_comb = pd.concat([df_rogfp2, df_oxd_rogfp2, df_rogfp2_cu], axis=0)
    df_comb.to_parquet("../data/sim-features.parquet", index=False)

    # Scale all features (excluding the system_label column)
    features = df_comb.drop(columns=["system_label"])
    scaler = MinMaxScaler()
    scaled_features = scaler.fit_transform(features)

    # Create a new DataFrame for the scaled data
    df_scaled = pd.DataFrame(scaled_features, columns=features.columns)
    df_scaled["system_label"] = df_comb[
        "system_label"
    ].values  # Ensure correct alignment
    df_scaled.to_parquet("../data/sim-scaled-features.parquet", index=False)

    # Store columns
    feat_cols = sorted(df_comb.columns.tolist())
    lines = [f"-   {col}\n" for col in feat_cols]
    with open("../sim-features.md", "w", encoding="utf-8") as f:
        f.writelines(lines)


if __name__ == "__main__":
    main()
