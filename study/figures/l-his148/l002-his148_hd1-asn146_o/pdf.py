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
        "analysis/005-rogfp-glh-md/data/struct-desc/his146_hd1-asn144_o-dist.npy",
    )
    rogfp_data = np.load(rogfp_data_path)
    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/his146_hd1-asn144_o-dist.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    # Copper
    rogfp_cu_data_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/his146_hd1-asn144_o-dist.npy",
    )
    rogfp_cu_data = np.load(rogfp_cu_data_path)

    # Compute all pdfs
    x_bounds = (1, 10)
    bin_width = 0.01  # Angstrom
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.05  # Manually tuned
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
    extrema_order = 5
    polyorder = 3
    window_length = int(0.5 / bin_width)
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Distance (Å)",
            pdf_rogfp_oxd,
            "Density",
            sci_notation=True,
            extrema_order=extrema_order,
            polyorder=polyorder,
            window_length=window_length,
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
    fig_title = "l002-his148_hd1-asn146_o"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "His148 HD1 - Asn146 O Distance [Å]"
    plot_x_bounds = (1, 6)
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
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_oxd, pmf_rogfp_cu = compute_pmfs(
        x_values, 3.22, (pdf_rogfp, pdf_rogfp_oxd, pdf_rogfp_cu), T=300.0
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
            extrema_order=extrema_order,
            polyorder=polyorder,
            window_length=window_length,
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
    plot_x_bounds = (1, 6)
    plot_y_bounds = (-1, 2)
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
