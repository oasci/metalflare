#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from metalflare.analysis.figures import use_mpl_rc_params
from metalflare.analysis.pdfs import extrema_table, make_pdf_fig

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
        "analysis/005-rogfp-glh-md/data/pocket/povme/volumes.csv",
    )
    rogfp_data = pd.read_csv(rogfp_data_path)["volume"].to_numpy()

    # Oxidized
    rogfp_oxd_data_path = os.path.join(
        base_dir,
        "analysis/007-rogfp-oxd-glh-md/data/pocket/povme/volumes.csv",
    )
    rogfp_oxd_data = pd.read_csv(rogfp_oxd_data_path)["volume"].to_numpy()

    rogfp2_cu_path = os.path.join(
        base_dir,
        "analysis/006-rogfp-cu-glh-md/data/pocket/povme/volumes.csv",
    )
    rogfp2_cu_data = pd.read_csv(rogfp2_cu_path)["volume"].to_numpy()

    rogfp2_na_path = os.path.join(
        base_dir,
        "analysis/008-rogfp-na-glh-md/data/pocket/povme/volumes.csv",
    )
    rogfp2_na_data = pd.read_csv(rogfp2_na_path)["volume"].to_numpy()

    x_bounds = (10, 500)
    x_values = np.arange(*x_bounds, 0.1)
    bw_method = 0.07

    kde = gaussian_kde(rogfp_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp_oxd_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_oxd = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp2_cu_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_cu = kde(x_values) / scaling_factor

    kde = gaussian_kde(rogfp2_na_data, bw_method=bw_method)
    scaling_factor = kde.integrate_box_1d(*x_bounds)
    pdf_rogfp_na = kde(x_values) / scaling_factor

    # save pdf information
    pdf_info_lines = ["Reduced roGFP2\n"]
    pdf_info_lines.extend(
        extrema_table(
            x_values,
            "Pocket volume",
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
            "Pocket volume",
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
            "Pocket volume",
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
            "Pocket volume",
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
    fig_title = "cro66_vol"
    pdf_plt_kwargs = {"alpha": 1.0, "linewidth": 1.5}
    x_label = "Cro66 Pocket Volume [Å$^3$]"
    plot_x_bounds = (75, 275)
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
    plt.xticks(np.arange(100, 250 + 0.01, 50))
    pdf_fig.savefig(f"{fig_title}-pdf.svg")
    plt.close()
