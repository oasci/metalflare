#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

os.chdir(os.path.dirname(os.path.realpath(__file__)))

rogfp2_file_path = "../data/rogfp2-desc.parquet"
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
            df.drop(col, axis=1, inplace=True)

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

    df_rogfp2_cu = pd.read_parquet(rogfp2_cu_file_path)
    df_rogfp2_cu = feat_angles(df_rogfp2_cu)

    # Scale all features
    scaler = MinMaxScaler()
    scaler.fit(np.vstack((df_rogfp2.values, df_rogfp2_cu.values)))

    rogfp2_arr = scaler.transform(df_rogfp2.values)
    np.save("../data/rogfp2-features.npy", rogfp2_arr)
    rogfp2_cu_arr = scaler.transform(df_rogfp2_cu.values)
    np.save("../data/rogfp2-cu-features.npy", rogfp2_cu_arr)


if __name__ == "__main__":
    main()
