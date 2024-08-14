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


def compute_pmfs(x, zero_at, pdfs, T=300.0):
    pmfs = []
    for pdf in pdfs:
        pmf = -KB * T * np.log(np.where(pdf > 0, pdf, 1e-50))
        pmf -= np.interp(zero_at, x, pmf)
        pmfs.append(pmf)
    return pmfs


def extrema_table(
    x,
    x_label,
    y,
    y_label,
    sci_notation=False,
    extrema_order=3,
    polyorder=2,
    window_length=4,
):
    local_extrema_x, local_extrema = get_extrema(
        x,
        y,
        extrema_order=extrema_order,
        polyorder=polyorder,
        window_length=window_length,
    )
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
    x_values,
    pdf_rogfp,
    pdf_rogfp_cu,
    plt_kwargs,
    x_label,
    x_bounds,
    y_label,
    y_bounds,
    pdf_rogfp_oxd=None,
    fill_between=False,
    legend_frame=False,
    figsize=(3.5, 3.0),
    pdf_na=None
):
    plt.figure(figsize=figsize)
    if fill_between:
        plot_function = plt.fill_between
    else:
        plot_function = plt.plot
    plot_function(
        x_values, pdf_rogfp, color="#1e2e79", label="Reduced", zorder=0, **plt_kwargs
    )
    if pdf_rogfp_oxd is not None:
        plot_function(
            x_values,
            pdf_rogfp_oxd,
            color="#EC4067",
            label="Oxidized",
            zorder=1,
            **plt_kwargs,
        )
    plot_function(
        x_values,
        pdf_rogfp_cu,
        color="#f99752",
        label="Cu(I)",
        zorder=2,
        **plt_kwargs,
    )
    if pdf_na is not None:
        plot_function(
            x_values,
            pdf_na,
            color="#1b998b",
            label="Na$^+$",
            alpha=0.6,
            linewidth=1.5,
            zorder=3
        )
    plt.xlabel(x_label)
    plt.xlim(*x_bounds)
    plt.ylim(*y_bounds)
    plt.ylabel(y_label)
    plt.legend(frameon=legend_frame)
    plt.tight_layout()
    return plt.gcf()


def make_pmf_fig(
    x_values,
    pmf_rogfp,
    pmf_rogfp_cu,
    x_label,
    x_bounds,
    y_label,
    y_bounds,
    pmf_rogfp_oxd=None,
    legend_frame=False,
    figsize=(3.5, 3.0)
):
    plt.figure(figsize=figsize)
    plt.plot(
        x_values, pmf_rogfp, color="#1e2e79", label="Reduced", zorder=0, linewidth=1.5
    )
    if pmf_rogfp_oxd is not None:
        plt.plot(
            x_values,
            pmf_rogfp_oxd,
            color="#EC4067",
            label="Oxidized",
            zorder=1,
            linewidth=1.5,
        )
    plt.plot(
        x_values,
        pmf_rogfp_cu,
        color="#f99752",
        label="Cu(I)",
        zorder=2,
        linewidth=1.5,
    )
    plt.axhline(y=0, linewidth=1.5, color="#C0C0C0", linestyle="dotted", zorder=-1)
    plt.xlabel(x_label)
    plt.xlim(*x_bounds)
    plt.ylim(*y_bounds)
    plt.ylabel(y_label)
    plt.legend(frameon=legend_frame)
    plt.tight_layout()
    return plt.gcf()
