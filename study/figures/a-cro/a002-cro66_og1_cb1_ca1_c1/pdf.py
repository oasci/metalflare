#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pdfs import (
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

    rogfp_data_path = os.path.join(
        base_dir,
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_og1_cb1_ca1_c1-dihedral.npy",
    )
    rogfp_data = np.load(rogfp_data_path)
    rogfp_data = np.degrees(rogfp_data)
    rogfp_data = np.concatenate([rogfp_data, rogfp_data + 360, rogfp_data - 360])

    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_og1_cb1_ca1_c1-dihedral.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_oxd_data = np.degrees(rogfp_oxd_data)
    rogfp_oxd_data = np.concatenate(
        [rogfp_oxd_data, rogfp_oxd_data + 360, rogfp_oxd_data - 360]
    )

    rogfp2_cu_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_og1_cb1_ca1_c1-dihedral.npy",
    )
    rogfp2_cu_data = np.load(rogfp2_cu_path)
    rogfp2_cu_data = np.degrees(rogfp2_cu_data)
    rogfp2_cu_data = np.concatenate(
        [rogfp2_cu_data, rogfp2_cu_data + 360, rogfp2_cu_data - 360]
    )

    x_bounds = (-180, 180)
    x_values = np.linspace(*x_bounds, 360 * 2)
    bw_method = 0.003

    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_oxd = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp2_cu_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_cu = kde(x_values) / scaling_factor

    # save pdf information
    pdf_info_lines = ["Reduced roGFP2\n"]
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Dihedral [°]",
            pdf_rogfp,
            "Density",
            sci_notation=True,
            window_length=20,
        )
    )
    pdf_info_lines.append("\nOxidized roGFP2\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Dihedral [°]",
            pdf_rogfp_oxd,
            "Density",
            sci_notation=True,
            window_length=20,
        )
    )
    pdf_info_lines.append("\nroGFP2 and Cu(I)\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Dihedral [°]",
            pdf_rogfp_cu,
            "Density",
            sci_notation=True,
            window_length=20,
        )
    )
    pdf_info_lines = [line + "\n" for line in pdf_info_lines]
    pdf_info_path = "./pdf-info.md"
    with open(pdf_info_path, "w", encoding="utf-8") as f:
        f.writelines(pdf_info_lines)

    # Make pdf plot
    fig_title = "a002-cro66_og1_cb1_ca1_c1"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 2.5}
    x_label = "Cro66 OG1-CB1-CA1-C1 Dihedral [°]"
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
        pdf_rogfp_oxd=pdf_rogfp_oxd,
    )
    plt.xticks(np.arange(-180, 181, 60))
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()

    # Compute potential of mean forces
    pmf_rogfp, pmf_rogfp_oxd, pmf_rogfp_cu = compute_pmfs(
        x_values, 48.32, (pdf_rogfp, pdf_rogfp_oxd, pdf_rogfp_cu), T=300.0
    )

    # save pmf information
    pmf_info_lines = ["Reduced roGFP2\n"]
    pmf_info_lines.extend(
        extrema_table(
            x_values, "Dihedral [°]", pmf_rogfp, "PMF [kcal/mol]", sci_notation=False
        )
    )
    pmf_info_lines.append("\nOxidized roGFP2\n")
    pmf_info_lines.extend(
        extrema_table(
            x_values,
            "Dihedral [°]",
            pmf_rogfp_oxd,
            "PMF [kcal/mol]",
            sci_notation=False,
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
    plot_y_bounds = (-4, 6)
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
    plt.xticks(np.arange(-180, 181, 60))
    pmf_fig.savefig(f"{fig_title}-pmf.svg")
    plt.close()
