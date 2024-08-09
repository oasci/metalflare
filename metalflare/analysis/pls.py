import string

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

from ..utils import format_feature_name


def perform_pls_regression(X, y, n_components=2):
    scaler_X = StandardScaler()
    scaler_y = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)
    y_scaled = y

    pls = PLSRegression(n_components=n_components, scale=False)
    pls.fit(X_scaled, y_scaled)

    return pls, X_scaled, y_scaled, scaler_X, scaler_y


def plot_pls_results(
    pls,
    X_scaled,
    y_scaled,
    feature_names,
    state_key,
    n_bins=150,
    y_label="Response variable",
    arrow_scaling=4.0,
    vmin=None,
    vmax=None,
):
    # Calculate PLS components
    pls_components = pls.transform(X_scaled)

    # Create the plot
    fig, ax = plt.subplots(figsize=(6, 5))

    # Create 2D histogram with average distance
    hist, xedges, yedges, im = ax.hist2d(
        pls_components[:, 0],
        pls_components[:, 1],
        bins=n_bins,
        weights=y_scaled,
        cmap="viridis",
    )

    # Normalize the histogram by the count in each bin to get the average
    counts, _, _ = np.histogram2d(
        pls_components[:, 0], pls_components[:, 1], bins=n_bins
    )
    hist_avg = np.divide(hist, counts, where=counts != 0)
    hist_avg[counts == 0] = np.nan  # Avoid division by zero

    # Update the image with the average values
    im.set_array(hist_avg.T)

    # Update colorbar to show average distance
    plt.colorbar(im, label=y_label)

    # Set the colorbar limits to the actual data range
    v_min, v_max = np.nanmin(hist_avg), np.nanmax(hist_avg)
    if vmin is not None:
        v_min = vmin
    if vmax is not None:
        v_max = vmax
    im.set_clim(v_min, v_max)

    # Add labels and title
    ax.set_xlabel("PLS Component 1")
    ax.set_ylabel("PLS Component 2")

    # Create a dictionary to map features to letters
    feature_letters = dict(
        zip(feature_names, string.ascii_uppercase[: len(feature_names)])
    )

    offset_distance = 0.5

    for label, loading in zip(feature_names, pls.x_loadings_):
        letter = feature_letters[label]
        loading *= arrow_scaling
        ax.arrow(
            0,
            0,
            loading[0],
            loading[1],
            color="#ff686b",
            alpha=1.0,
            head_width=0.11,
            head_length=0.1,
            width=0.03,
        )

        # Calculate unit vector in direction of arrow
        magnitude = np.sqrt(loading[0] ** 2 + loading[1] ** 2)
        unit_vector = loading / magnitude if magnitude != 0 else np.array([0, 0])

        # Calculate label position
        label_x = loading[0] + unit_vector[0] * offset_distance
        label_y = loading[1] + unit_vector[1] * offset_distance

        ax.annotate(
            letter,
            (label_x, label_y),
            xytext=(0, 0),
            textcoords="offset points",
            ha="center",
            va="center",
            alpha=1.0,
            fontweight="bold",
            bbox=dict(boxstyle="circle,pad=0.3", fc="white", ec="none", alpha=0.75),
        )

    # Add R-squared value
    r_squared = pls.score(X_scaled, y_scaled)
    ax.text(
        0.05,
        0.95,
        f"R² = {r_squared:.3f}",
        transform=ax.transAxes,
        verticalalignment="top",
    )

    # Calculate the gradient of y values
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(x_centers, y_centers)

    # Use np.gradient to compute the gradient of hist_avg
    dy, dx = np.gradient(hist_avg.T)

    # Calculate the average gradient, ignoring NaN values
    avg_dx = np.nanmean(dx)
    avg_dy = np.nanmean(dy)

    gradient_vector = None

    if np.isfinite(avg_dx) and np.isfinite(avg_dy):
        # Normalize the gradient vector
        gradient_magnitude = np.sqrt(avg_dx**2 + avg_dy**2)
        if gradient_magnitude > 0:
            normalized_dx = avg_dx / gradient_magnitude
            normalized_dy = avg_dy / gradient_magnitude
            gradient_vector = np.array([normalized_dx, normalized_dy])

            # Calculate the endpoints of the line
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            # Find the intersections of the line with the plot borders
            t_values = [
                (xlim[0] - 0) / normalized_dx,  # Left border
                (xlim[1] - 0) / normalized_dx,  # Right border
                (ylim[0] - 0) / normalized_dy,  # Bottom border
                (ylim[1] - 0) / normalized_dy,  # Top border
            ]

            # Sort t_values to get the endpoints
            t_min, t_max = min(t_values), max(t_values)

            # Calculate the endpoints
            x1, y1 = t_min * normalized_dx, t_min * normalized_dy
            x2, y2 = t_max * normalized_dx, t_max * normalized_dy

            # Plot the dashed line
            ax.plot([x1, x2], [y1, y2], color="#AFAFAF", linestyle="--", linewidth=2)

        else:
            print("Warning: Gradient magnitude is zero. Skipping gradient line.")
    else:
        print("Warning: Invalid gradient values. Skipping gradient line.")

    plt.tight_layout()
    plt.savefig(f"{state_key}_pls_regression.png", dpi=300)
    plt.close()

    return gradient_vector


def compute_loading_angles_and_magnitudes(pls, gradient_vector, feature_names):
    loadings = pls.x_loadings_[:, :2]  # Only consider the first two components
    magnitudes = np.linalg.norm(loadings, axis=1)
    scaler = StandardScaler()
    magnitudes = scaler.fit_transform(np.array(magnitudes).reshape(-1, 1)).flatten()

    sin_angles = []
    for loading in loadings:
        # Compute the sine of the angle between the loading and the gradient vector
        cos_angle = np.dot(loading, gradient_vector) / (
            np.linalg.norm(loading) * np.linalg.norm(gradient_vector)
        )
        sin_angle = np.sqrt(1 - cos_angle**2)  # sin^2 + cos^2 = 1
        # Determine the sign of the sine
        cross_product = np.cross(loading, gradient_vector)
        sin_angle = np.sign(cross_product) * sin_angle
        sin_angles.append(sin_angle)

    # Create a list of tuples (feature, sine of angle, magnitude)
    result = list(zip(feature_names, sin_angles, magnitudes))

    result.sort(key=lambda x: -abs(x[2]))

    return result


def create_angle_magnitude_file(data, state_key, feature_names):
    feature_letters = dict(
        zip(feature_names, string.ascii_uppercase[: len(feature_names)])
    )

    output_file = f"{state_key}_loading_angles_and_magnitudes.md"

    with open(output_file, "w") as f:
        f.write(r"| Feature | $\sin \left( \theta \right)$ | Magnitude |")
        f.write("\n|---------|------------|-----------|\n")

        for feature, sin_angle, magnitude in data:
            letter = feature_letters[feature]
            formatted_feature = format_feature_name(feature)
            f.write(
                f"| **{letter}**: {formatted_feature} | {sin_angle:.4f} | {magnitude:.4f} |\n"
            )


def create_markdown_table(data, headers, feature_names):
    """Create a markdown table from the given data and headers, including feature letters."""
    # Create a dictionary to map features to letters
    feature_letters = dict(
        zip(feature_names, string.ascii_uppercase[: len(feature_names)])
    )

    table = "| " + " | ".join(headers) + " |\n"
    table += "|" + "---|" * len(headers) + "\n"
    for row in data:
        feature = row[0]
        letter = feature_letters[feature]
        formatted_feature = format_feature_name(feature)
        formatted_row = [f"**{letter}**: {formatted_feature}"] + row[1:]
        table += "| " + " | ".join(map(str, formatted_row)) + " |\n"
    return table


def compare_states(results, feature_names, output_file="loadings_analysis.md"):
    # Extract loadings and create a DataFrame
    loadings_df = {}
    for state in results:
        loadings = results[state][0].x_loadings_
        # Calculate magnitude using the first two components
        magnitudes = np.sqrt(loadings[:, 0] ** 2 + loadings[:, 1] ** 2)
        scaler = StandardScaler()
        magnitudes = scaler.fit_transform(np.array(magnitudes).reshape(-1, 1)).flatten()
        loadings_df[state] = magnitudes
    loadings_df = pd.DataFrame(loadings_df)
    loadings_df.index = feature_names

    sort_by = loadings_df.mean(axis=1)
    sorted_df = loadings_df.loc[sort_by.sort_values(ascending=False).index]

    # Prepare data for the table
    data = []
    for feature, row in sorted_df.iterrows():
        reduced = row["reduced"]
        oxidized = row["oxidized"]
        cu = row["cu"]
        oxidized_diff = oxidized - reduced
        cu_diff = cu - reduced

        data.append(
            [
                feature,
                f"{reduced:.4f}",
                f"{oxidized:.4f} ({oxidized_diff:+.4f})",
                f"{cu:.4f} ({cu_diff:+.4f})",
            ]
        )

    # Write to file
    with open(output_file, "w") as f:
        f.write(
            create_markdown_table(
                data, ["Feature", "Reduced", "Oxidized", "Cu"], feature_names
            )
        )
        f.write("\n")
