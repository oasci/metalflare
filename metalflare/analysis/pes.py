import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

KB = 1.987204259e-3  # kcal/(mol K)
"""Boltzmann constant in kcal/(mol K)."""


def get_2d_hist(
    x_data: ArrayLike,
    y_data: ArrayLike,
    bins: int | ArrayLike | tuple[int, int] | tuple[ArrayLike, ArrayLike] = 13,
    T: float | int | None = None,
) -> tuple[np.ma.MaskedArray, dict[str, tuple[ArrayLike, ArrayLike]]]:
    r"""Compute a 2D histogram and calculate the potential of mean force.

    This function creates a 2D histogram from the input data and calculates
    the potential of mean force (PMF) or the negative logarithm of probability,
    depending on whether a temperature is provided.

    Args:
        x_data: 1D array of x-coordinates for the data points.
        y_data: 1D array of y-coordinates for the data points.
            Must have the same length as x_data.
        bins: The bin specification for the 2D histogram.
            - If int, the number of bins for both dimensions.
            - If array-like, the bin edges for both dimensions.
            - If [int, int], the number of bins in each dimension.
            - If [array, array], the bin edges in each dimension.
        T: Temperature in Kelvin. If provided, the function
            calculates the potential of mean force. If `None`, it calculates
            $-ln(p)$ where $p$ is the probability.

    Returns:
        2D array containing the calculated histogram values. Zero counts are masked.
            If `T` is provided, values represent the potential of mean force in energy
            units $k_B T$. If `T` is `None`, values represent $-\ln(p)$.
        A dictionary containing information about the
            bin edges and centers:

            - 'edges': tuple of `(x_edges, y_edges)`
            - 'centers': tuple of `(x_centers, y_centers)`

    Raises:
        ValueError: If `x_data` and `y_data` have different lengths.

    Notes:
        - The function uses numpy.histogram2d to create the initial 2D histogram.
        - Zero counts in the histogram are masked to avoid log(0) errors.
        - The histogram is normalized to represent probabilities.

    Example:
        ```python
        >>> import numpy as np
        >>> x = np.random.randn(1000)
        >>> y = np.random.randn(1000)
        >>> hist, bins_info = get_2d_hist(x, y, bins=20, T=300)
        >>> print(hist.shape)
        (20, 20)
        >>> print(bins_info.keys())
        dict_keys(['edges', 'centers'])
        ```
    """
    hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=bins)
    hist = np.ma.masked_where(hist == 0, hist)
    hist /= np.sum(hist)
    hist = -np.log(hist)
    if T is not None:
        hist *= KB * T
        hist -= np.min(hist)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    bins_info = {"edges": (x_edges, y_edges), "centers": (x_centers, y_centers)}
    return hist, bins_info


def create_pes(
    x_data: ArrayLike,
    y_data: ArrayLike,
    bins: int | ArrayLike | tuple[int, int] | tuple[ArrayLike, ArrayLike] = 30,
    vmin: float = 0,
    vmax: float = 10,
    levels: int = 30,
    T: float | None = None,
    ax: mpl.axis.Axis | None = None,
    colorbar: bool = True,
) -> plt.Figure:
    """Create a potential energy surface (PES) plot using matplotlib.

    This function generates a 2D histogram from the input data using the get_2d_hist function,
    and then creates a contour plot representing either the potential of mean force (PMF)
    or the negative logarithm of probability.

    Args:
        x_data: 1D array of x-coordinates for the data points.
        y_data: 1D array of y-coordinates for the data points.
        bins: The bin specification for the 2D histogram. Defaults to 30.

            - If int, the number of bins for both dimensions.
            - If array-like, the bin edges for both dimensions.
            - If (int, int), the number of bins in each dimension.
            - If (array, array), the bin edges in each dimension.
        vmin: Minimum value for the colorbar. Defaults to 0.
        vmax: Maximum value for the colorbar. Defaults to 10.
        levels: Number of contour levels. Defaults to 30.
        T: Temperature in Kelvin. If provided, the plot represents PMF.
            If None, it represents -ln(p). Defaults to None.

    Returns:
        A matplotlib Figure object containing the PES plot.

    Raises:
        ValueError: If x_data and y_data have different lengths.

    Example:
        ```python
        >>> import numpy as np
        >>> x = np.random.randn(1000)
        >>> y = np.random.randn(1000)
        >>> fig = create_pes(x, y, bins=20, T=300)
        >>> fig.savefig('pes_plot.png')
        ```
    """
    hist, bins_info = get_2d_hist(x_data, y_data, bins=bins, T=T)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    if ax is None:
        ax = plt.gca()  # Get the current axis if no axis is provided

    contour = ax.contourf(
        *bins_info["centers"], hist.T, levels=levels, cmap="viridis", norm=norm
    )

    if T is not None:
        label = "PMF [kcal/mol]"
    else:
        label = "-ln(p)"

    if colorbar:
        ax.figure.colorbar(
            contour,
            ax=ax,
            label=label,
            norm=norm,
            ticks=list(range(vmin, int(np.floor(vmax)) + 1)),
        )

    ax.set_aspect("auto")
    return ax.figure


def masked_difference(
    hist_ref: np.ma.MaskedArray, hist: np.ma.MaskedArray
) -> np.ma.MaskedArray:
    """
    Compute the difference between two masked histograms with special handling for masked values.

    This function calculates the difference between two histograms, taking into account
    masked values. It handles special cases where data is present in one histogram but
    not in the other.

    Args:
        hist_ref: The reference histogram.
        hist: The histogram to compare against the reference.

    Returns:
        The difference between the two histograms, with special
            handling for masked values:

            -   Where both histograms have values, it returns the normal
                `hist` - `hist_ref` difference.
            -   Where `hist_ref` is masked and `hist` is not, it returns the difference
                between `hist`'s value and the maximum value in `hist`.
            -   Where `hist` is masked and `hist_ref` is not, it returns `hist_ref`'s
                value.

    Raises:
        ValueError: If the input arrays have different shapes.

    Note:
        This function assumes that masked values in the histograms represent
        configurations that were not observed. The special handling of these cases
        is designed to capture significant changes in the observability of states.
    """
    diff = hist - hist_ref

    # Handle the special cases
    mask_not_in_ref = np.where(hist_ref.mask & ~hist.mask)
    mask_only_in_ref = np.where(~hist_ref.mask & hist.mask)
    diff[mask_not_in_ref] = hist[mask_not_in_ref] - np.max(hist)
    diff[mask_only_in_ref] = hist_ref[mask_only_in_ref]

    return diff


def create_pes_difference(
    data_x_ref: ArrayLike,
    data_y_ref: ArrayLike,
    data_x: ArrayLike,
    data_y: ArrayLike,
    bins: int | ArrayLike | tuple[int, int] | tuple[ArrayLike, ArrayLike] = 30,
    vmin: float = -10,
    vmax: float = 10,
    levels: int = 30,
    T: float | None = None,
) -> plt.Figure:
    """Create a difference plot of two potential energy surfaces (PES) using matplotlib.

    This function generates 2D histograms from two sets of input data using the get_2d_hist function,
    calculates their difference, and then creates a contour plot representing either the difference
    in potential of mean force (ΔPMF) or the difference in negative logarithm of probability.

    Args:
        data_x_ref: 1D array of x-coordinates for the reference data points.
        data_y_ref: 1D array of y-coordinates for the reference data points.
        data_x: 1D array of x-coordinates for the comparison data points.
        data_y: 1D array of y-coordinates for the comparison data points.
        bins: The bin specification for the 2D histograms.

            - If int, the number of bins for both dimensions.
            - If array-like, the bin edges for both dimensions.
            - If (int, int), the number of bins in each dimension.
            - If (array, array), the bin edges in each dimension.
        vmin: Minimum value for the colorbar.
        vmax: Maximum value for the colorbar.
        levels: Number of contour levels.
        T: Temperature in Kelvin. If provided, the plot represents ΔPMF.
            If None, it represents -Δln(p).

    Returns:
        A matplotlib Figure object containing the PES difference plot.

    Raises:
        ValueError: If input data arrays have different lengths or are incompatible.

    Example:
        ```python
        >>> import numpy as np
        >>> x_ref = np.random.randn(1000)
        >>> y_ref = np.random.randn(1000)
        >>> x = np.random.randn(1000)
        >>> y = np.random.randn(1000)
        >>> fig = create_pes_difference(x_ref, y_ref, x, y, bins=20, T=300)
        >>> fig.savefig('pes_difference_plot.png')
        ```

    Notes:
        - The function uses a TwoSlopeNorm for coloring, centered at 0.
        - The resulting plot uses a diverging colormap (RdBu_r) to highlight differences.
        - The masked_difference function should be defined elsewhere to handle
          the difference calculation between potentially masked arrays.
    """
    hist_ref, bins_info = get_2d_hist(data_x_ref, data_y_ref, bins=bins, T=T)
    hist, _ = get_2d_hist(data_x, data_y, bins=bins_info["edges"], T=T)

    hist_diff = masked_difference(hist_ref, hist)

    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    contour = plt.contourf(
        *bins_info["centers"], hist_diff.T, levels=levels, cmap="RdBu_r", norm=norm
    )

    if T is not None:
        label = "ΔPMF [kcal/mol]"
    else:
        label = "-Δln(p)"
    plt.colorbar(
        contour,
        label=label,
        norm=norm,
        ticks=list(range(vmin, int(np.floor(vmax)) + 1)),
    )

    plt.tight_layout()
    return plt.gcf()
