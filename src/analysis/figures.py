"""This module contains the plotting functions for the article."""

import glob

import matplotlib.pyplot as plt
import numpy as np

from src.analysis.cvs import strings_to_SF_IG
from src.analysis.plotting import (add_XRD_values, plot_2D_heatmap,
                                   two_cv_strings_time_series)
from src.analysis.utils import natural_sort

plt.rcParams["axes.facecolor"] = "#f9f9fb"
plt.rcParams["grid.color"] = "white"
plt.rcParams["grid.linestyle"] = "-"
plt.rcParams["grid.linewidth"] = 2
plt.rcParams["axes.grid"] = True
plt.rcParams["lines.solid_capstyle"] = "round"


def final_FES_IG_SF(
    name,
    path_processed,
    path_report,
    XRD_dictionary,
    fig_title,
    show_cbar=False,
    version="",
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
    fig.savefig(f"{path_report}/FES_{name}{version}.png")
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
    fig.savefig(f"{path_report}/FES_{cv_name}_path.png")
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
    fig.savefig(f"{path_report}/projection_{cv_name}.png")
    return fig, ax


def final_2D_string_convergence(
    name,
    path_report,
    fig_title,
):
    files = natural_sort(glob.glob(f"{name}/strings/string[0-9]*txt"))
    strings = np.array([np.loadtxt(file).T for file in files])
    reduced_string = strings_to_SF_IG(strings, [0, 1], [10, 11])
    reduced_string_labels = ["SF (nm)", "IG (nm)"]
    fig, ax = two_cv_strings_time_series(
        reduced_string,
        reduced_string_labels,
        start_iteration=0,
        n_average=50,
        av_last_n_it=50,
        fig_title=fig_title,
        position="lower left",
    )
    fig.tight_layout()
    fig.savefig(f"{path_report}/convergence_2D_{name}.png")
    return fig, ax


def final_1D_path(path_processed, name, error, label, color, path_report, version=""):
    fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    F = np.load(f"{path_processed}{name}/FES_path.npy")
    ax.plot(F[0], F[1], marker="", label=label, color=color, lw=3)
    n_boot = error[0]
    error_block = error[1]
    e = np.load(f"{path_processed}{name}/path_errors_{n_boot}_5.npy")
    e = e[error_block, :]
    ax.fill_between(F[0], F[1] + e, F[1] - e, alpha=0.1, color=color)
    ax.legend(prop={"size": 18}, loc="best")
    ax.set_xlabel(r"$s_{path}$", size=25)
    ax.set_ylabel("Free Energy ($k_BT$)", size=20)
    ax.set_ylim([-0.5, 20])
    ax.set_xlim([0.1, 1.01])
    ax.tick_params("x", labelsize=15)
    ax.tick_params("y", labelsize=15)
    ax.grid(False)
    fig.savefig(path_report + f"FES_path_cv_{name}{version}.png")

    return fig, ax
