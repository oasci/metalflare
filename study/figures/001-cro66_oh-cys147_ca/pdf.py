#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

os.chdir(os.path.dirname(os.path.realpath(__file__)))


def find_local_maxima_1d(arr):
    """
    Finds the local maxima in a 1D numpy array.

    Parameters:
    - arr: A 1D numpy array

    Returns:
    - A tuple containing two elements:
        - A list of local maxima values
        - A list of positions (indices) of the local maxima
    """
    # Ensure the array is 1D
    if arr.ndim != 1:
        raise ValueError("The function expects a 1D numpy array.")

    # Shift the array to the right and left
    right_shifted = np.roll(arr, 1)
    left_shifted = np.roll(arr, -1)

    # Identify local maxima
    local_max_mask = (arr > right_shifted) & (arr > left_shifted)

    # Handle boundary conditions by setting the first and last element to False
    local_max_mask[0] = local_max_mask[-1] = False

    # Extract local maxima values and their positions
    local_maxima = arr[local_max_mask]
    positions = np.where(local_max_mask)[0]

    sort_idx = np.argsort(local_maxima)[::-1]
    local_maxima = local_maxima[sort_idx]
    positions = positions[sort_idx]

    return local_maxima, positions


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../"

    rogfp_dist_path = os.path.join(
        base_dir, "analysis/001-rogfp-md/data/struct-desc/cro65_oh-cym145_ca-dist.npy"
    )
    rogfp_dist = np.load(rogfp_dist_path)

    rogfp_cu_dist_path = os.path.join(
        base_dir,
        "analysis/003-rogfp-cu-md/data/struct-desc/cro65_oh-cym145_ca-dist.npy",
    )
    rogfp_cu_dist = np.load(rogfp_cu_dist_path)

    kwargs = {"alpha": 0.5, "linewidth": 1.0}
    x_bounds = (1, 10)
    x_values = np.linspace(*x_bounds, 1000)
    bw_method = 0.1

    print("\nrogfp2")
    kde = gaussian_kde(rogfp_dist, bw_method=bw_method)
    rogfp_kde_values = kde(x_values)
    local_maxima, positions = find_local_maxima_1d(rogfp_kde_values)
    print("y-axis local maxima:", local_maxima)
    print("x-values:", x_values[positions])
    plt.fill_between(
        x_values, rogfp_kde_values, color="#1e2e79", label="Unbound", **kwargs
    )

    print("\nrogfp2-cu")
    kde = gaussian_kde(rogfp_cu_dist, bw_method=bw_method)
    rogfp_cu_kde_values = kde(x_values)
    local_maxima, positions = find_local_maxima_1d(rogfp_cu_kde_values)
    print("y-axis local maxima:", local_maxima)
    print("x-values:", x_values[positions])
    plt.fill_between(
        x_values, rogfp_cu_kde_values, color="#f99752", label="Bound", **kwargs
    )

    plt.xlabel("CRO66 OH - CYS147 CA Distance [Å]")
    plt.xlim(4.5, 9)
    plt.ylim(0, None)
    plt.ylabel("Density")

    plt.legend()
    plt.tight_layout()

    plt.savefig("001-cro66_oh-cys147_ca-pdf.svg")
    plt.close()

    # Use boltzmann inversion
    kb = 1.987204259e-3  # kcal/(mol K)
    T = 300.0  # K

    pmf_rogfp = (
        -kb * T * np.log(np.where(rogfp_kde_values > 0, rogfp_kde_values, 1e-15))
    )
    pmf_rogfp_cu = (
        -kb * T * np.log(np.where(rogfp_cu_kde_values > 0, rogfp_cu_kde_values, 1e-15))
    )
    shift_by = min(np.min(pmf_rogfp), np.min(pmf_rogfp_cu))
    pmf_rogfp -= shift_by
    pmf_rogfp_cu -= shift_by

    plt.plot(x_values, pmf_rogfp, color="#1e2e79", label="Unbound", linewidth=2.5)
    plt.plot(x_values, pmf_rogfp_cu, color="#f99752", label="Bound", linewidth=2.5)
    plt.xlabel("CRO66 OH - CYS147 CA Distance [Å]")
    plt.xlim(4.5, 9)
    plt.ylim(0, 5)
    plt.ylabel("PMF [kcal/mol]")
    plt.xlabel("CRO66 OH - THR203 OG1 Distance [Å]")
    plt.legend()
    plt.tight_layout()

    plt.savefig("001-cro66_oh-cys147_ca-pmf.svg")
    plt.close()
