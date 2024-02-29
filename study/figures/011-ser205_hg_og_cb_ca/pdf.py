#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.pdfs import (
    compute_pmfs,
    extrema_table,
    make_pdf_fig,
    make_pmf_fig,
)

os.chdir(os.path.dirname(os.path.realpath(__file__)))


if __name__ == "__main__":
    # Specify the paths to the trajectory and topology files
    base_dir = "../../"

    rogfp_data_path = os.path.join(
        base_dir,
        "analysis/001-rogfp-md/data/struct-desc/ser203_hg_og_cb_ca-dihedral.npy",
    )
    rogfp_data = np.load(rogfp_data_path)
    rogfp_data = np.degrees(rogfp_data)
    rogfp_data = np.concatenate([rogfp_data, rogfp_data + 360, rogfp_data - 360])

    rogfp2_cu_path = os.path.join(
        base_dir,
        "analysis/003-rogfp-cu-md/data/struct-desc/ser203_hg_og_cb_ca-dihedral.npy",
    )
    rogfp2_cu_data = np.load(rogfp2_cu_path)
    rogfp2_cu_data = np.degrees(rogfp2_cu_data)
    rogfp2_cu_data = np.concatenate(
        [rogfp2_cu_data, rogfp2_cu_data + 360, rogfp2_cu_data - 360]
    )

    x_bounds = (-180, 180)
    x_values = np.linspace(*x_bounds, 1000)
    bw_method = 0.03

    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(-180, 180)
    pdf_rogfp = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp2_cu_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(-180, 180)
    pdf_rogfp_cu = kde(x_values) / scaling_factor

    # save pdf information
    pdf_info_lines = ["roGFP2\n"]
    pdf_info_lines.extend(
        extrema_table(x_values, "Dihedral [°]", pdf_rogfp, "Density", sci_notation=True)
    )
    pdf_info_lines.append("\nroGFP2 and Cu(I)\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values, "Dihedral [°]", pdf_rogfp_cu, "Density", sci_notation=True
        )
    )
    pdf_info_lines = [line + "\n" for line in pdf_info_lines]
    pdf_info_path = "./pdf-info.md"
    with open(pdf_info_path, "w", encoding="utf-8") as f:
        f.writelines(pdf_info_lines)

    # Make pdf plot
    fig_title = "011-ser205_hg_og_cb_ca"
    pdf_plt_kwargs = {"alpha": 0.5, "linewidth": 1.0}
    x_label = "SER205 HG-OH-CB-CA Dihedral [°]"
    plot_x_bounds = (-180, 180)
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
    )
    plt.xticks(np.arange(-180, 181, 60))
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_cu = compute_pmfs(
        pdf_rogfp, pdf_rogfp_cu, x_values, -66.13, T=300.0
    )

    # save pmf information
    pmf_info_lines = ["roGFP2\n"]
    pmf_info_lines.extend(
        extrema_table(
            x_values, "Dihedral [°]", pmf_rogfp, "PMF [kcal/mol]", sci_notation=False
        )
    )
    pmf_info_lines.append("\nroGFP2 and Cu(I)\n")
    pmf_info_lines.extend(
        extrema_table(
            x_values, "Dihedral [°]", pmf_rogfp_cu, "PMF [kcal/mol]", sci_notation=False
        )
    )
    pmf_info_lines = [line + "\n" for line in pmf_info_lines]
    pmf_info_path = "./pmf-info.md"
    with open(pmf_info_path, "w", encoding="utf-8") as f:
        f.writelines(pmf_info_lines)

    y_label = "PMF [kcal/mol]"
    plot_y_bounds = (None, None)
    pmf_fig = make_pmf_fig(
        x_values,
        pmf_rogfp,
        pmf_rogfp_cu,
        x_label=x_label,
        x_bounds=plot_x_bounds,
        y_label=y_label,
        y_bounds=plot_y_bounds,
    )
    plt.xticks(np.arange(-180, 181, 60))
    pmf_fig.savefig(f"{fig_title}-pmf.svg")
    plt.close()
