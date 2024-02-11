"""Assist with making probability distribution functions."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from ..utils import get_extrema

KB = 1.987204259e-3  # kcal/(mol K)


def compute_pdf(arr1, x_values, bw_method=None):
    kde = gaussian_kde(arr1, bw_method=bw_method)
    pdf = kde(x_values)
    return pdf


def compute_pmfs(pdf1, pdf2, T=300.0):
    pmf1 = -KB * T * np.log(np.where(pdf1 > 0, pdf1, 1e-15))
    pmf2 = -KB * T * np.log(np.where(pdf2 > 0, pdf2, 1e-15))
    shift_by = min(np.min(pmf1), np.min(pmf2))
    pmf1 -= shift_by
    pmf2 -= shift_by
    return pmf1, pmf2


def extrema_table(x, x_label, y, y_label, sci_notation=False):
    local_extrema_x, local_extrema = get_extrema(x, y)
    lines = []
    # Start of the Markdown table
    lines.append(f"| {x_label} | {y_label} |")
    lines.append("|-----------|-----------|")
    # Loop through the values to format and print each row
    for x_val, y_val in zip(local_extrema_x, local_extrema):
        if sci_notation:
            y_val_str = f"${y_val:.3e}".replace("e", " \\times 10^{") + "}$"
            y_val_str = y_val_str.replace("e+0", "").replace("e-0", "-")
        else:
            y_val_str = f"{y_val:.3f}"
        lines.append(f"| {x_val:.2f} | {y_val_str} |")
    return lines


def make_pdf_fig(
    x_values, pdf_rogfp, pdf_rogfp_cu, plt_kwargs, x_label, x_bounds, y_label, y_bounds
):
    plt.fill_between(
        x_values, pdf_rogfp, color="#1e2e79", label="Unbound", **plt_kwargs
    )
    plt.fill_between(
        x_values, pdf_rogfp_cu, color="#f99752", label="Bound", **plt_kwargs
    )
    plt.xlabel(x_label)
    plt.xlim(*x_bounds)
    plt.ylim(*y_bounds)
    plt.ylabel(y_label)
    plt.legend()
    plt.tight_layout()
    return plt.gcf()


def make_pmf_fig(
    x_values, pmf_rogfp, pmf_rogfp_cu, x_label, x_bounds, y_label, y_bounds
):
    plt.plot(x_values, pmf_rogfp, color="#1e2e79", label="Unbound", linewidth=2.5)
    plt.plot(x_values, pmf_rogfp_cu, color="#f99752", label="Bound", linewidth=2.5)
    plt.xlabel(x_label)
    plt.xlim(*x_bounds)
    plt.ylim(*y_bounds)
    plt.ylabel(y_label)
    plt.legend()
    plt.tight_layout()
    return plt.gcf()
