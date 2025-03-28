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
        "analysis/005-rogfp-glh-md/data/struct-desc/cro65_n3_ca3_c3-val66_n-dihedral.npy",
    )
    rogfp_data = np.load(rogfp_data_path)
    rogfp_data = np.degrees(rogfp_data)
    rogfp_data = np.concatenate([rogfp_data, rogfp_data + 360, rogfp_data - 360])

    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/struct-desc/cro65_n3_ca3_c3-val66_n-dihedral.npy",
    )
    rogfp_oxd_data = np.load(rogfp_oxd_data_path)
    rogfp_oxd_data = np.degrees(rogfp_oxd_data)
    rogfp_oxd_data = np.concatenate(
        [rogfp_oxd_data, rogfp_oxd_data + 360, rogfp_oxd_data - 360]
    )

    rogfp2_cu_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/struct-desc/cro65_n3_ca3_c3-val66_n-dihedral.npy",
    )
    rogfp2_cu_data = np.load(rogfp2_cu_path)
    rogfp2_cu_data = np.degrees(rogfp2_cu_data)
    rogfp2_cu_data = np.concatenate(
        [rogfp2_cu_data, rogfp2_cu_data + 360, rogfp2_cu_data - 360]
    )

    rogfp2_na_path = os.path.join(
        base_dir,
        "analysis/008-rogfp-na-glh-md/data/struct-desc/cro65_n3_ca3_c3-val66_n-dihedral.npy",
    )
    rogfp2_na_data = np.load(rogfp2_na_path)
    rogfp2_na_data = np.degrees(rogfp2_na_data)
    rogfp2_na_data = np.concatenate(
        [rogfp2_na_data, rogfp2_na_data + 360, rogfp2_na_data - 360]
    )

    x_bounds = (-120, 240)
    x_values = np.linspace(*x_bounds, 360 * 2)
    bw_method = 0.003

    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp = kde(x_values) / scaling_factor
    reduced_fraction = kde.integrate_box_1d(120, 240) / scaling_factor
    print(f"Reduced kde stat:  {reduced_fraction:.3f}")

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_oxd = kde(x_values) / scaling_factor
    oxidized_fraction = kde.integrate_box_1d(120, 240) / scaling_factor
    print(f"Oxidized kde stat:  {oxidized_fraction:.3f}")

    kde = gaussian_kde(rogfp2_cu_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_cu = kde(x_values) / scaling_factor
    cu_fraction = kde.integrate_box_1d(120, 240) / scaling_factor
    print(f"Cu(I) kde stat:  {cu_fraction:.3f}")

    kde = gaussian_kde(rogfp2_na_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_na = kde(x_values) / scaling_factor
    na_fraction = kde.integrate_box_1d(120, 240) / scaling_factor
    print(f"Na+ kde stat:  {na_fraction:.3f}")

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
    pdf_info_lines.append("\nroGFP2 and Na+\n")
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Dihedral [°]",
            pdf_rogfp_na,
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
    fig_title = "a004-cro66_n3_ca3_c3_val68_n"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "Cro66 $\psi$ [deg.]"
    plot_x_bounds = (-120, 240)
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
        pdf_na=pdf_rogfp_na,
    )
    plt.xticks(np.arange(-120, 240 + 0.1, 60))
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()
