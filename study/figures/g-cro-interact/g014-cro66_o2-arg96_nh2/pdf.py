#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
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
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_o2-arg94_nh2-dist.npy",
    )
    rogfp_data = np.load(rogfp_dist_path)
    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_o2-arg94_nh2-dist.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_cu_dist_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_o2-arg94_nh2-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_dist_path)

    rogfp_na_dist_path = os.path.join(
        base_dir,
        "analysis/008-rogfp-na-glh-md/data/struct-desc/cro65_o2-arg94_nh2-dist.npy",
    )
    rogfp_na_data = np.load(rogfp_na_dist_path)

    # Compute all pdfs
    x_bounds = (1, 10)
    bin_width = 0.05  # Angstrom
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.01  # Manually adjusted
    pdf_rogfp = compute_pdf(rogfp_data, x_values, bw_method=bw_method)
    pdf_rogfp_oxd = compute_pdf(rogfp_oxd_data, x_values, bw_method=bw_method)
    pdf_rogfp_cu = compute_pdf(rogfp_cu_data, x_values, bw_method=bw_method)
    pdf_rogfp_na = compute_pdf(rogfp_na_data, x_values, bw_method=bw_method)

    # KDE stats
    x_min = 1
    x_max = 4
    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    reduced_fraction = kde.integrate_box_1d(x_min, x_max)
    print(f"Reduced kde stat:  {reduced_fraction:.3f}")

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    oxidized_fraction = kde.integrate_box_1d(x_min, x_max)
    print(f"Oxidized kde stat: {oxidized_fraction:.3f}")

    kde = gaussian_kde(rogfp_cu_data, bw_method=bw_method)
    cu_fraction = kde.integrate_box_1d(x_min, x_max)
    print(f"Cu(I) kde stat:    {cu_fraction:.3f}")

    kde = gaussian_kde(rogfp_na_data, bw_method=bw_method)
    na_fraction = kde.integrate_box_1d(x_min, x_max)
    print(f"Na+ kde stat:     {na_fraction:.3f}")

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
    fig_title = "g014-cro66_o2-arg96_nh2"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "Cro66 O2 - Arg95 NH2 Distance [Å]"
    plot_x_bounds = (2, 4)
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
        figsize=(3.5, 3.0),
        pdf_na=pdf_rogfp_na,
    )

    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()
