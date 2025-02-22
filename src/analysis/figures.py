"""This module contains the plotting functions for the article."""

import glob
import pickle

import matplotlib.pyplot as plt
import numpy as np

from src.analysis.cvs import strings_to_SF_IG
from src.analysis.plotting import (
    add_XRD_values,
    all_rmsd_strings_time_series,
    plot_2D_heatmap,
    plot_trajectories_map,
    strings_time_series,
    two_cv_strings_time_series,
)
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
    XRD_dictionary=None,
    fig_title=None,
    show_cbar=False,
    version="",
    restarts=False,
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
    if XRD_dictionary is not None:
        add_XRD_values(
            XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left"
        )
    fig.tight_layout()
    ax.set_box_aspect(1)
    if restarts:
        trajectories = np.load(f"{path_processed}/{name}/SF_IG_restarts.npy")
        ax = plot_trajectories_map(ax, trajectories)

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
    version="",
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
    ax.set_box_aspect(1)
    fig.savefig(f"{path_report}/FES_{cv_name}_path{version}.png")
    return fig, ax


def final_cv_projection(
    name,
    path_processed,
    path_report,
    cv_name,
    cv_data,
    cv_label,
    XRD_dictionary=None,
    fig_title=None,
    f_max=None,
    f_min=None,
    show_cbar=False,
    fig=None,
    ax=None,
    cmap=plt.cm.viridis_r,
):
    F = np.load(f"{path_processed}/{name}/FES_SF_IG.npy")
    extent = np.load(f"{path_processed}/{name}/extent.npy")

    if fig is None:
        fig, ax = plot_2D_heatmap(
            cv_data,
            extent,
            cmap=cmap,
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
        add_XRD_values(
            XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left"
        )
        ax.set_xlim([0.48, 1.01])
        ax.set_ylim([1.1, 2.45])
        ax.set_box_aspect(1)
        fig.tight_layout()
        fig.savefig(f"{path_report}/projection_{cv_name}.png")
    else:
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
            fig=fig,
            ax=ax,
        )
        add_XRD_values(
            XRD_dictionary, "SF", "IG", size=15, ax=ax, position="lower left"
        )
        ax.set_xlim([0.48, 1.01])
        ax.set_ylim([1.1, 2.45])
        ax.set_box_aspect(1)
    return fig, ax


def final_2D_string_convergence(
    name, path_data, path_report, fig_title, inset_pos=None
):
    files = natural_sort(glob.glob(f"{path_data}/{name}/strings/string[0-9]*txt"))
    strings = np.array([np.loadtxt(file).T for file in files])
    reduced_string = strings_to_SF_IG(strings, [0, 1], [10, 11])
    reduced_string_labels = ["SF (nm)", "IG (nm)"]
    fig, aa = two_cv_strings_time_series(
        reduced_string,
        reduced_string_labels,
        start_iteration=0,
        n_average=50,
        av_last_n_it=50,
        fig_title=fig_title,
        position="lower left",
    )
    ax2 = fig.add_axes(inset_pos)
    two_cv_strings_time_series(
        reduced_string[:10, :, :],
        [None, None],
        start_iteration=0,
        n_average=1,
        av_last_n_it=1,
        fig_title=None,
        position="lower left",
        ax=ax2,
        fig=fig,
        inset=True,
    )
    aa.set_box_aspect(1)
    ax2.set_box_aspect(1)
    fig.tight_layout()
    fig.savefig(f"{path_report}/convergence_2D_{name}.png")


def final_1D_path(
    path_processed, name, error, label, color, path_report, version="", annotations=None
):
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
    ax.set_box_aspect(1)
    if annotations is not None:
        for position, n in zip(annotations["positions"], annotations["names"]):
            ax.text(
                position[0],
                position[1],
                n,
                fontsize=annotations["fontsize"],
                rotation=0,
                color=annotations["color"],
                va="center",
                ha="center",
            )
    fig.savefig(path_report + f"FES_path_cv_{name}{version}.png")

    return fig, ax


def final_rmsd_string(path_data, path_report, sim_name, version=""):
    fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    colors = ["C0", "C1", "C3"]
    for name, c in zip(sim_name.keys(), colors):
        files = natural_sort(glob.glob(f"{path_data}/{name}/strings/string[0-9]*txt"))
        strings = np.array([np.loadtxt(file).T for file in files])
        label = sim_name[name]
        fig, ax = all_rmsd_strings_time_series(
            strings, "RMSD[String] (nm)", label=label, fig=fig, ax=ax, color=c
        )
    ax.tick_params("x", labelsize=15)
    ax.tick_params("y", labelsize=15)
    ax.grid(False)
    fig.savefig(path_report + f"RMSD_string{version}.png")
    return fig, ax


def final_FES_IG_SF_error(
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
    errors = np.load(f"{path_processed}/{name}/errors_150_5.npy")
    e = errors[3]
    e[~np.isfinite(F)] = np.nan

    e_max = 4
    f_max = 20
    fig, ax = plot_2D_heatmap(
        e,
        extent,
        cmap=plt.cm.viridis_r,
        f_max=e_max,
        f_min=0,
        c_density=F,
        c_max=f_max,
        cbar_label="Error ($k_BT$)",
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
    ax.set_box_aspect(1)
    fig.savefig(f"{path_report}/FES_error_{name}{version}.png")
    return fig, ax


def final_SF_content(name, path_processed, path_report, fig_title, version=""):
    fig, ax = plt.subplots(4, 2, figsize=(10, 5 * 4), sharex=True, sharey=True)
    color_maps = (plt.cm.Blues, plt.cm.Purples)
    colors = ("k", "k")
    F = np.load(f"{path_processed}/{name}/FES_SF_IG.npy")
    extent = np.load(f"{path_processed}/{name}/extent.npy")
    with open(f"{path_processed}/{name}/SF_occupation_data.pkl", "rb") as handle:
        SF_occupation_data = pickle.load(handle)
    for i, ligand in enumerate(("W", "K")):
        for j in range(4):
            _ = plot_2D_heatmap(
                SF_occupation_data[ligand][j + 1],
                extent,
                # cbar_label=f"Number of {ligand} in S{j}",
                cbar_label=f"Number of {ligand}",
                xlabel="SF (nm)",
                ylabel="IG (nm)",
                f_min=0,
                f_max=1,
                ax=ax[j, i],
                fig=fig,
                cmap=color_maps[i],
                n_colors=200,
                c_density=F,
                c_min=0,
                c_max=25,
                c_color=colors[i],
                labelsize=13,
            )
            ax[j, i].grid(None)
            ax[j, i].set_title(f"Site {j+1}", pad=15, size=15)
    fig.tight_layout()
    fig.savefig(f"{path_report}/SF_content_per_site_{name}{version}.png")
    return fig, ax


def final_strings_vs_time_series(
    name, path_data, path_report, fig_title, select, ylabels, version=""
):
    with open(f"{path_data}/{name}/cv.pkl", "rb") as file:
        cvs, ndx_groups = pickle.load(file)
    ndx_groups = {x: ndx_groups[x] for x in select}
    files = natural_sort(glob.glob(f"{path_data}/{name}/strings/string[0-9]*txt"))
    strings = np.array([np.loadtxt(file).T for file in files])
    strings = strings[:, -4:, :]
    fig, ax = strings_time_series(
        strings,
        ndx_groups,
        start_iteration=1,
        n_average=25,
        av_last_n_it=None,
        sharex=True,
        sharey=True,
    )
    for i, a in enumerate(ax.flatten()):
        a.set_ylabel(ylabels[i])
        a.set_title("")
    fig.tight_layout()
    fig.savefig(f"{path_report}/L81-W67_string_time_series_{name}{version}.png")

    return fig, ax


def final_cv_vs_time_series(
    name,
    path_data,
    system_name,
    path_report,
    fig_title,
    select,
    ylabels,
    starts,
    ends,
    version="",
):
    n_beads = 16
    strings = np.load(f"{path_data}/{system_name}/{name}.npy")[::64, :]
    it = int(strings.shape[0] / n_beads)
    strings = np.reshape(strings, [it, n_beads, 4], order="C")
    strings = np.hstack([strings, ends * np.ones([it, 1, 4])])
    strings = np.hstack([starts * np.ones([it, 1, 4]), strings])
    strings = np.swapaxes(strings, 1, 2)
    fig, ax = strings_time_series(
        strings,
        None,
        start_iteration=1,
        n_average=25,
        av_last_n_it=None,
        sharex=True,
        sharey=True,
        y_axis_tic=0.2,
    )
    for i, a in enumerate(ax.flatten()):
        a.set_ylabel(ylabels[i])
        a.set_title("")
    fig.tight_layout()
    fig.savefig(f"{path_report}/{name}_string_time_series_{system_name}{version}.png")

    return fig, ax
