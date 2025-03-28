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

    # Reduced
    rogfp_data_path = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/thr201_og1-glu220_he2-dist.npy",
    )
    rogfp_data = np.load(rogfp_data_path)

    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/thr201_og1-glu220_he2-dist.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)

    # Copper
    rogfp_cu_data_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/thr201_og1-glu220_he2-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_data_path)

    rogfp_na_dist_path = os.path.join(
        base_dir,
        "analysis/008-rogfp-na-glh-md/data/struct-desc/thr201_og1-glu220_he2-dist.npy",
    )
    rogfp_na_data = np.load(rogfp_na_dist_path)

    # Compute all pdfs
    x_bounds = (1, 15)
    bin_width = 0.1  # Angstrom
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.05  # Manually tuned
    pdf_rogfp = compute_pdf(rogfp_data, x_values, bw_method=bw_method)
    pdf_rogfp_oxd = compute_pdf(rogfp_oxd_data, x_values, bw_method=bw_method)
    pdf_rogfp_cu = compute_pdf(rogfp_cu_data, x_values, bw_method=bw_method)
    pdf_rogfp_na = compute_pdf(rogfp_na_data, x_values, bw_method=bw_method)

    # KDE stats
    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    reduced_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Reduced kde stat:  {reduced_fraction:.3e}")

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    oxidized_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Oxidized kde stat: {oxidized_fraction:.3e}")

    kde = gaussian_kde(rogfp_cu_data, bw_method=bw_method)
    cu_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Cu(I) kde stat:    {cu_fraction:.3e}")

    kde = gaussian_kde(rogfp_na_data, bw_method=bw_method)
    na_fraction = kde.integrate_box_1d(0.1, 2.5)
    print(f"Na+ kde stat:      {na_fraction:.3e}")

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
    fig_title = "d007-thr203_og1-glu222_he2"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "Thr203 OG1 - Glu222 HE2 Distance [Å]"
    plot_x_bounds = (1, 10)
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

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_oxd, pmf_rogfp_cu, pmf_rogfp_na = compute_pmfs(
        x_values, 7.65, (pdf_rogfp, pdf_rogfp_oxd, pdf_rogfp_cu, pdf_rogfp_na), T=300.0
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
    plot_x_bounds = (1, 10)
    plot_y_bounds = (-1, 6)
    pmf_fig = make_pmf_fig(
        x_values,
        pmf_rogfp,
        pmf_rogfp_cu,
        x_label=x_label,
        x_bounds=plot_x_bounds,
        y_label=y_label,
        y_bounds=plot_y_bounds,
        pmf_rogfp_oxd=pmf_rogfp_oxd,
        pmf_rogfp_na=pmf_rogfp_na,
    )
    pmf_fig.savefig(f"{fig_title}-pmf.svg")
    plt.close()
