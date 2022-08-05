"""This module contains the plotting functions for the article."""

import matplotlib.pyplot as plt
import numpy as np
from src.analysis.plotting import plot_2D_heatmap, add_XRD_values


plt.rcParams["axes.facecolor"] = "#f9f9fb"
plt.rcParams["grid.color"] = "white"
plt.rcParams["grid.linestyle"] = "-"
plt.rcParams["grid.linewidth"] = 2
plt.rcParams["axes.grid"] = True
plt.rcParams["lines.solid_capstyle"] = "round"


def final_FES_IG_SF(
    name, path_processed, path_report, XRD_dictionary, fig_title, show_cbar=False
):

    F = np.load(f"{path_processed}/{name}/FES_SF_IG.npy")
    extent = np.load(f"{path_processed}/{name}/extent.npy")

    f_max = 20
    fig, ax = plot_2D_heatmap(
        F,
        extent,
        f_max=f_max,
        f_min=0,
        cbar_label="Free Energy ($k_BT$)",
        xlabel="Selectivity Filter (nm)",
        ylabel="Inner Gate (nm)",
        fig_title=fig_title,
        show_grid=False,
        show_cbar=show_cbar,
    )
    ax.set_xlim([0.48, 1.01])
    ax.set_ylim([1.1, 2.45])
    add_XRD_values(XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left")
    fig.tight_layout()
    fig.savefig(f"{path_report}/main_panel/FES_{name}.png")
    return fig, ax


def final_FES_path_CV(
    name,
    cv_name,
    path_processed,
    path_report,
    XRD_dictionary,
    fig_title,
    ylabel,
    show_cbar=False,
):

    F = np.load(f"{path_processed}/{name}/FES_{cv_name}_path.npy")
    extent = np.load(f"{path_processed}/{name}/extent_{cv_name}_path.npy")

    f_max = 20
    fig, ax = plot_2D_heatmap(
        F,
        extent,
        f_max=f_max,
        f_min=0,
        cbar_label="Free Energy ($k_BT$)",
        xlabel=r"$s_{path}$",
        ylabel=ylabel,
        fig_title=fig_title,
        show_grid=False,
        show_cbar=show_cbar,
    )
    # ax.set_xlim([0.48, 1.01])
    # ax.set_ylim([1.1, 2.45])
    # add_XRD_values(XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left")
    fig.tight_layout()
    fig.savefig(f"{path_report}/path_vs_CV/FES_{cv_name}_path.png")
    return fig, ax


def final_cv_projection(
    name,
    path_processed,
    path_report,
    cv_name,
    cv_data,
    cv_label,
    XRD_dictionary,
    fig_title,
    f_max=None,
    f_min=None,
    show_cbar=False,
):

    F = np.load(f"{path_processed}/{name}/FES_SF_IG.npy")
    extent = np.load(f"{path_processed}/{name}/extent.npy")

    fig, ax = plot_2D_heatmap(
        cv_data,
        extent,
        cmap=plt.cm.viridis_r,
        f_max=f_max,
        f_min=f_min,
        cbar_label=cv_label,
        xlabel="Selectivity Filter (nm)",
        ylabel="Inner Gate (nm)",
        fig_title=fig_title,
        show_grid=False,
        show_cbar=show_cbar,
        c_density=F,
    )
    add_XRD_values(XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left")
    ax.set_xlim([0.48, 1.01])
    ax.set_ylim([1.1, 2.45])
    fig.tight_layout()
    fig.savefig(f"{path_report}/path_vs_CV/projection_{cv_name}.png")
    return fig, ax
