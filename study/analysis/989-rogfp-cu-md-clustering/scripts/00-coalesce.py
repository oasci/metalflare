#!/usr/bin/env python3

import glob
import os

import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.realpath(__file__)))

rogfp2_dset_path = "../../001-rogfp-md/data/struct-desc"
rogfp2_cu_dset_path = "../../003-rogfp-cu-md/data/struct-desc"


def npy_to_parquet(file_paths, output_file):
    """
    Reads multiple .npy files, converts them into a Pandas DataFrame, and saves it as a Parquet file.

    :param file_paths: List of file paths to .npy files.
    :param output_file: Path to the output Parquet file.
    """
    df_list = (
        []
    )  # Initialize an empty list to store DataFrames created from .npy files.
    df_columns = []

    for file_path in file_paths:
        file_name = os.path.basename(file_path.replace(".npy", ""))

        # Load the numpy array from .npy file.
        np_array = np.load(file_path)

        df_columns.append(file_name)

        # Convert numpy array to Pandas DataFrame.
        df = pd.DataFrame(np_array)

        # Append the DataFrame to our list.
        df_list.append(df)

    # Concatenate all DataFrames in the list into a single DataFrame.
    final_df = pd.concat(df_list, ignore_index=True, axis=1)
    final_df.columns = df_columns

    # Save the DataFrame as a Parquet file.
    final_df.to_parquet(output_file, index=False)


def main():
    rogfp2_file_paths = glob.glob(os.path.join(rogfp2_dset_path, "*.npy"))
    file_path = "../data/rogfp2-desc.parquet"
    npy_to_parquet(rogfp2_file_paths, file_path)

    rogfp2_cu_file_paths = glob.glob(os.path.join(rogfp2_cu_dset_path, "*.npy"))
    file_path = "../data/rogfp2-cu-desc.parquet"
    npy_to_parquet(rogfp2_cu_file_paths, file_path)


if __name__ == "__main__":
    main()
