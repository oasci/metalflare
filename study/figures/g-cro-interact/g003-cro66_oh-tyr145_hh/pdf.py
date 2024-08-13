#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pdfs import (
    compute_pdf,
    compute_pmfs,
    extrema_table,
    make_pdf_fig,
    make_pmf_fig,
)

os.chdir(os.path.dirname(os.path.realpath(__file__)))


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../../"

    # Update plot params
    rc_json_path = os.path.join(
        base_dir, "misc/003-figure-style/matplotlib-rc-params.json"
    )
    font_dirs = [os.path.join(base_dir, "misc/003-figure-style/roboto")]
    use_mpl_rc_params(rc_json_path, font_dirs)

    rogfp_dist_path = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_oh-tyr143_hh-dist.npy",
    )
    rogfp_data = np.load(rogfp_dist_path)
    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_oh-tyr143_hh-dist.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_cu_dist_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_oh-tyr143_hh-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_dist_path)

    # Compute all pdfs
    x_bounds = (1, 10)
    bin_width = 0.05  # Angstrom
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.04  # Manually adjusted
    pdf_rogfp = compute_pdf(rogfp_data, x_values, bw_method=bw_method)
    pdf_rogfp_oxd = compute_pdf(rogfp_oxd_data, x_values, bw_method=bw_method)
    pdf_rogfp_cu = compute_pdf(rogfp_cu_data, x_values, bw_method=bw_method)

    # KDE stats
    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    reduced_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Reduced kde stat:  {reduced_fraction:.3f}")

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    oxidized_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Oxidized kde stat: {oxidized_fraction:.3f}")

    kde = gaussian_kde(rogfp_cu_data, bw_method=bw_method)
    cu_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Cu(I) kde stat:    {cu_fraction:.3f}")

    # save pdf information
    pdf_info_lines = ["Reduced roGFP2\n"]
    pdf_info_lines.extend(
        extrema_table(x_values, "Distance (Å)", pdf_rogfp, "Density", sci_notation=True)
    )
    pdf_info_lines.append("\nOxidized roGFP2\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values, "Distance (Å)", pdf_rogfp_oxd, "Density", sci_notation=True
        )
    )
    pdf_info_lines.append("\nroGFP2 and Cu(I)\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values, "Distance (Å)", pdf_rogfp_cu, "Density", sci_notation=True
        )
    )
    pdf_info_lines = [line + "\n" for line in pdf_info_lines]
    pdf_info_path = "./pdf-info.md"
    with open(pdf_info_path, "w", encoding="utf-8") as f:
        f.writelines(pdf_info_lines)

    # Make pdf plot
    fig_title = "g003-cro66_oh-tyr145_hh"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "Cro66 OH - Tyr145 HH Distance [Å]"
    plot_x_bounds = (1, 7)
    y_label = "Density"
    plot_y_bounds = (0, None)

    pdf_fig = make_pdf_fig(
        x_values,
        pdf_rogfp,
        pdf_rogfp_cu,
        pdf_plt_kwargs,
        x_label=x_label,
        x_bounds=plot_x_bounds,
        y_label=y_label,
        y_bounds=plot_y_bounds,
        pdf_rogfp_oxd=pdf_rogfp_oxd,
    )

    # Add hydrogen bond region
    rect = Rectangle((0, 0), 2.15, 10, facecolor="#F5F5F5", zorder=-10)
    plt.gca().add_patch(rect)
    colors = ["#F5F5F5", "#ffffff"]
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    plt.imshow(gradient, extent=[2.14, 2.5, 0, 10], aspect="auto", cmap=cmap, zorder=-9)

    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_oxd, pmf_rogfp_cu = compute_pmfs(
        x_values, 3.21, (pdf_rogfp, pdf_rogfp_oxd, pdf_rogfp_cu), T=300.0
    )

    # save pmf information
    pmf_info_lines = ["Reduced roGFP2\n"]
    pmf_info_lines.extend(
        extrema_table(
            x_values, "Distance (Å)", pmf_rogfp, "PMF [kcal/mol]", sci_notation=False
        )
    )
    pmf_info_lines.append("\nOxidized roGFP2\n")
    pmf_info_lines.extend(
        extrema_table(
            x_values,
            "Distance (Å)",
            pmf_rogfp_oxd,
            "PMF [kcal/mol]",
            sci_notation=False,
        )
    )
    pmf_info_lines.append("\nroGFP2 and Cu(I)\n")
    pmf_info_lines.extend(
        extrema_table(
            x_values, "Distance (Å)", pmf_rogfp_cu, "PMF [kcal/mol]", sci_notation=False
        )
    )
    pmf_info_lines = [line + "\n" for line in pmf_info_lines]
    pmf_info_path = "./pmf-info.md"
    with open(pmf_info_path, "w", encoding="utf-8") as f:
        f.writelines(pmf_info_lines)

    y_label = "PMF [kcal/mol]"
    plot_y_bounds = (-2, 2)
    pmf_fig = make_pmf_fig(
        x_values,
        pmf_rogfp,
        pmf_rogfp_cu,
        x_label=x_label,
        x_bounds=plot_x_bounds,
        y_label=y_label,
        y_bounds=plot_y_bounds,
        pmf_rogfp_oxd=pmf_rogfp_oxd,
    )
    pmf_fig.savefig(f"{fig_title}-pmf.svg")
    plt.close()
