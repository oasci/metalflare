#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np

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

    # A distance of 0.0 means no water was found.
    # Need to replace 0.0 with higher number.

    # Reduced
    rogfp_data_path_1 = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_oh-h2o_h1-dist.npy",
    )
    rogfp_data_path_2 = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_oh-h2o_h2-dist.npy",
    )
    rogfp_data_1 = np.load(rogfp_data_path_1)
    rogfp_data_1[rogfp_data_1 == 0.0] = 20.0
    rogfp_data_2 = np.load(rogfp_data_path_2)
    rogfp_data_2[rogfp_data_2 == 0.0] = 20.0
    rogfp_data = np.minimum(rogfp_data_1, rogfp_data_2)

    # Oxidized
    rogfp_oxd_data_path_1 = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_oh-h2o_h1-dist.npy",
    )
    rogfp_oxd_data_path_2 = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_oh-h2o_h2-dist.npy",
    )
    rogfp_oxd_data_1 = np.load(rogfp_oxd_data_path_1)
    rogfp_oxd_data_1[rogfp_oxd_data_1 == 0.0] = 20.0
    rogfp_oxd_data_2 = np.load(rogfp_oxd_data_path_2)
    rogfp_oxd_data_2[rogfp_oxd_data_2 == 0.0] = 20.0
    rogfp_oxd_data = np.minimum(rogfp_oxd_data_1, rogfp_oxd_data_2)

    # Copper
    rogfp_cu_data_path_1 = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_oh-h2o_h1-dist.npy",
    )
    rogfp_cu_data_path_2 = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_oh-h2o_h2-dist.npy",
    )
    rogfp_cu_data_1 = np.load(rogfp_cu_data_path_1)
    rogfp_cu_data_1[rogfp_cu_data_1 == 0.0] = 20.0
    rogfp_cu_data_2 = np.load(rogfp_cu_data_path_2)
    rogfp_cu_data_2[rogfp_cu_data_2 == 0.0] = 20.0
    rogfp_cu_data = np.minimum(rogfp_cu_data_1, rogfp_cu_data_2)

    # Compute all pdfs
    x_bounds = (0, 10)
    bin_width = 0.1  # Angstrom
    n_bins = int((max(x_bounds) - min(x_bounds)) / bin_width)
    x_values = np.linspace(*x_bounds, n_bins)
    bw_method = 0.005
    pdf_rogfp = compute_pdf(rogfp_data, x_values, bw_method=bw_method)
    pdf_rogfp_oxd = compute_pdf(rogfp_oxd_data, x_values, bw_method=bw_method)
    pdf_rogfp_cu = compute_pdf(rogfp_cu_data, x_values, bw_method=bw_method)

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
    fig_title = "e009-cro66_oh-h2o_h"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 2.5}
    x_label = "CRO66 OH - H2O H Distance [Å]"
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
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_oxd, pmf_rogfp_cu = compute_pmfs(
        x_values, 3.17, (pdf_rogfp, pdf_rogfp_oxd, pdf_rogfp_cu), T=300.0
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
    plot_x_bounds = (1, 7)
    plot_y_bounds = (-3, 2)
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
