"""This module contains the plotting functions and defaults of the project."""

from math import ceil

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

from src.analysis.string_tica_msm import get_kde

plt.rcParams["axes.facecolor"] = "#f9f9fb"
plt.rcParams["grid.color"] = "white"
plt.rcParams["grid.linestyle"] = "-"
plt.rcParams["grid.linewidth"] = 2
plt.rcParams["axes.grid"] = True
plt.rcParams["lines.solid_capstyle"] = "round"


def _colorbar(mappable, cmap, norm, label0, size=10):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    ax = mappable.axes
    # fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cbar.set_label(label0, size=size)
    return cbar


def plot_2D_heatmap(
    G,
    extent,
    cmap=None,
    f_min=None,
    f_max=None,
    cbar_label="",
    ax=None,
    fig=None,
    xlabel="",
    n_colors=50,
    ylabel="",
    c_density=None,
    levels=20,
    c_min=None,
    c_max=None,
    c_color="k",
    fig_title=None,
    show_grid=True,
    show_cbar=True,
):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    if cmap is None:
        cmap = plt.cm.RdYlBu_r
    if f_max is None:
        f_max = np.nanmax(G)
    if f_min is None:
        f_min = np.nanmin(G)
    if c_max is None:
        c_max = np.nanmax(c_density)
    if c_min is None:
        c_min = np.nanmin(c_density)
    ax.contourf(G, n_colors, cmap=cmap, extent=extent, vmin=f_min, vmax=f_max)
    norm = mpl.colors.Normalize(vmin=f_min, vmax=f_max)
    if show_cbar:
        _ = _colorbar(ax, cmap, norm, cbar_label, 20)
    if c_density is not None:
        ax.contour(
            c_density,
            levels=levels,
            extent=extent,
            vmin=c_min,
            vmax=c_max,
            colors=c_color,
        )
    ax.set_xlabel(xlabel, size=23, labelpad=15)
    ax.set_ylabel(ylabel, size=23, labelpad=15)
    ax.set_title(fig_title, size=27, pad=15)
    ax.tick_params("x", labelsize=15)
    ax.tick_params("y", labelsize=15)
    ax.grid(show_grid)
    fig.tight_layout()
    return fig, ax


def strings_time_series(
    strings, ndx_groups, start_iteration=0, n_average=1, av_last_n_it=None
):

    n_plots = strings.shape[1]
    n_strings = strings.shape[0]
    fig, ax = plt.subplots(ceil(n_plots / 2), 2, figsize=(20, 8 * ceil(n_plots / 2)))
    ax = ax.flatten()
    cmap = plt.cm.viridis_r
    n_colors = (n_strings - start_iteration) // n_average + 1
    colors = cmap(np.linspace(0, 1, n_colors))  # yellow to blue
    norm = mpl.colors.Normalize(vmin=start_iteration, vmax=n_strings - 1)

    for i, a in enumerate(ax[:n_plots]):
        a.plot(strings[0, i, :], ls=":", marker=".", label="string0", color="r")
        for jj, j in enumerate(range(start_iteration, n_strings, n_average)):
            string = np.mean(strings[j : j + n_average, i, :], axis=0)
            a.plot(string, ls="-", marker="o", color=colors[jj])
            av = np.mean(strings[n_strings - av_last_n_it :, i, :], axis=0)
        if av_last_n_it is not None:
            std = np.std(strings[n_strings - av_last_n_it :, i, :], axis=0)
            a.fill_between(
                np.arange(len(av)),
                av + std,
                av - std,
                alpha=0.4,
                color="mediumorchid",
                label=f"std(string{start_iteration}-{n_strings})",
            )

            a.plot(
                av,
                ls="-",
                lw=3,
                marker=".",
                color="mediumorchid",
                label=f"mean(string{start_iteration}-{n_strings})",
            )
        a.set_ylabel(
            f"{list(ndx_groups.keys())[2*i]} - {list(ndx_groups.keys())[2*i+1]} (nm)",
            size=18,
            labelpad=16,
        )
        a.set_xlabel("bead number", size=15, labelpad=13)
        a.set_xlim(left=0, right=strings.shape[2] - 1)
        a.xaxis.set_minor_locator(MultipleLocator(1))
        a.xaxis.set_major_locator(MultipleLocator(1))
        a.yaxis.set_minor_locator(MultipleLocator(0.1))
        a.yaxis.set_major_locator(MultipleLocator(0.1))
        a.grid(which="minor")
        a.tick_params(axis="y", labelsize=14)
        a.tick_params(axis="x", labelsize=11)
        a.set_title(f"cv{i}")
        if i % 2 != 0:
            a.legend()
            _colorbar(a, cmap, norm, "iteration number", 20)
    if n_plots % 2:
        fig.delaxes(ax[-1])
    return fig, ax


def rmsd_strings_time_series(strings, ndx_groups):
    n_plots = strings.shape[1]
    n_strings = strings.shape[0]
    fig, ax = plt.subplots(ceil(n_plots / 2), 2, figsize=(20, 8 * ceil(n_plots / 2)))
    ax = ax.flatten()
    for i, a in enumerate(ax[:n_plots]):
        x = np.arange(n_strings)
        y = strings[:, i, :] - strings[0, i, :]
        y = np.sqrt(np.sum(y * y, axis=1) / strings.shape[2])
        a.plot(x, y)

        a.set_ylabel(
            f"RMSD[{list(ndx_groups.keys())[2*i]} - {list(ndx_groups.keys())[2*i+1]} (nm)]",
            size=18,
            labelpad=16,
        )
        a.set_xlabel("iteration number", size=15, labelpad=13)
        a.set_title(f"cv{i}")
    if n_plots % 2:
        fig.delaxes(ax[-1])
    return fig, ax


def two_cv_strings_time_series(
    strings,
    strings_labels,
    start_iteration=0,
    n_average=1,
    av_last_n_it=None,
    fig_title=None,
    position="best",
):

    n_strings = strings.shape[0]
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    cmap = plt.cm.viridis_r
    n_colors = (n_strings - start_iteration) // n_average + 1
    colors = cmap(np.linspace(0, 1, n_colors))  # yellow to blue
    norm = mpl.colors.Normalize(vmin=start_iteration, vmax=n_strings - 1)
    ax.plot(
        strings[0, 0, :],
        strings[0, 1, :],
        ls=":",
        marker=".",
        label="string0",
        color="k",
    )
    if av_last_n_it is not None:
        av_0 = np.mean(strings[n_strings - av_last_n_it :, 0, :], axis=0)
        std_0 = np.std(strings[n_strings - av_last_n_it :, 0, :], axis=0)
        av_1 = np.mean(strings[n_strings - av_last_n_it :, 1, :], axis=0)
        std_1 = np.std(strings[n_strings - av_last_n_it :, 1, :], axis=0)
        ax.plot(
            av_0,
            av_1,
            ls="-",
            lw=3,
            marker=".",
            color="mediumorchid",
            label=f"mean(string{n_strings - av_last_n_it}-{n_strings})",
        )

    for jj, j in enumerate(range(start_iteration, n_strings, n_average)):
        av_0 = np.mean(strings[j:, 0, :], axis=0)
        std_0 = np.std(strings[j:, 0, :], axis=0)
        av_1 = np.mean(strings[j:, 1, :], axis=0)
        std_1 = np.std(strings[j:, 1, :], axis=0)
        ax.errorbar(
            av_0, av_1, fmt="--", xerr=std_0, yerr=std_1, color=colors[jj], alpha=0.9
        )
    ax.set_ylabel(
        strings_labels[1],
        size=23,
        labelpad=15,
    )
    ax.set_xlabel(
        strings_labels[0],
        size=23,
        labelpad=15,
    )

    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.grid(which="minor")
    ax.tick_params("x", labelsize=15)
    ax.tick_params("y", labelsize=15)
    ax.set_title(fig_title, size=27, pad=15)
    ax.legend(loc=position, prop={"size": 15})
    _colorbar(ax, cmap, norm, "iteration number", 20)
    return fig, ax


def all_rmsd_strings_time_series(strings, ylabel, label=None, fig=None, ax=None):
    n_strings = strings.shape[0]
    if fig is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    x = np.arange(n_strings)
    y = strings[:, :, :] - strings[0, :, :]
    y = np.sqrt(np.sum(y * y, axis=(1, 2)) / strings.shape[2])
    ax.plot(x, y, label=label)
    ax.set_ylabel(
        ylabel,
        size=18,
        labelpad=16,
    )
    ax.legend(loc="best", prop={"size": 15})
    ax.set_xlabel("iteration number", size=15, labelpad=13)
    return fig, ax


def plot_FES_1d_vs_t(FES_vs_t, xlabel=None, cmap=None, fig=None, ax=None, error=None):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    if cmap is None:
        cmap = plt.cm.viridis_r
    n_fes = FES_vs_t.shape[0] - 1
    colors = cmap(np.linspace(0, 1, n_fes))  # yellow to blue
    norm = mpl.colors.Normalize(vmin=0, vmax=n_fes)

    # ax.plot(FES_vs_t[0, :], FES_vs_t[-1, :], marker="o", c="red", ms=6, label="full")
    for i in range(1, n_fes + 1):
        ax.plot(FES_vs_t[0, :], FES_vs_t[i, :], marker=".", c=colors[i - 1])
    _ = _colorbar(ax, cmap, norm, "Iteration", 15)
    if error is not None:
        ax.fill_between(
            FES_vs_t[0, :], FES_vs_t[-1, :] + error, FES_vs_t[-1, :] - error, alpha=0.3
        )

    ax.set_ylabel("Free Energy (kT)", size=18)
    ax.set_xlabel(xlabel, size=18)
    ax.legend()

    return fig, ax


def add_XRD_values(
    XRD_dictionary,
    valx=None,
    valy=None,
    size=15,
    color="k",
    ax=None,
    fig=None,
    txt_size=None,
    linestyle="-",
    position="best",
):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    for key in XRD_dictionary.keys():
        state = XRD_dictionary[key]["state"]
        label = f"{key} ({state})"
        if "color" in XRD_dictionary[key].keys():
            color0 = XRD_dictionary[key]["color"]
        else:
            color0 = color
        if valx is not None and valy is not None:
            x = XRD_dictionary[key][valx]
            y = XRD_dictionary[key][valy]
            marker = XRD_dictionary[key]["marker"]
            ax.scatter(
                x,
                y,
                marker=marker,
                s=size**2,
                c=color0,
                label=label,
                # linewidth=2,
                # edgecolors="w",
            )
        elif valy is not None:
            y = XRD_dictionary[key][valy]
            if "linestyle" in XRD_dictionary[key].keys():
                linestyle = XRD_dictionary[key]["linestyle"]
            ax.axhline(
                y,
                linestyle=linestyle,
                c=color0,
                label=label,
                lw=size,
            )
        elif valx is not None:
            x = XRD_dictionary[key][valx]
            if "linestyle" in XRD_dictionary[key].keys():
                linestyle = XRD_dictionary[key]["linestyle"]
            ax.axvline(
                x,
                linestyle=linestyle,
                c=color0,
                label=label,
                lw=size,
            )
    if txt_size is None:
        txt_size = size
    ax.legend(loc=position, prop={"size": txt_size})

    return


def plot_path_vs_cv_FES(
    s_path,
    other_cv,
    weights,
    cv_name,
    cv_fig_label,
    bandwidth,
    f_max,
    save_npy,
    save_fig,
    save_extent,
):
    cvs = np.concatenate([s_path, other_cv], axis=2)
    bandwidth = 0.05
    p_of_cv, extent = get_kde(cvs, weights, bandwidth)
    F0 = -np.log(p_of_cv)
    F = F0 - F0.min()
    F[F > f_max] = np.nan
    np.save(save_npy, F)
    np.save(save_extent, np.array(extent))
    fig, ax = plot_2D_heatmap(
        F,
        extent,
        f_max=f_max,
        f_min=0,
        cbar_label="Free Energy (kT)",
        xlabel="s[path]",
        ylabel=cv_fig_label,
    )
    fig.tight_layout()
    fig.savefig(save_fig)
    return fig, ax


def plot_trajectories_map(ax, trajectories):
    for i, trajectory in enumerate(trajectories):
        ax.plot(
            trajectory[0] / 10,
            trajectory[1] / 10,
            ls="",
            marker="o",
            c=f"C{i}",
            ms=8,
            label=f"Simulation {i}",
        )
        ax.plot(
            trajectory[0][:1] / 10,
            trajectory[1][:1] / 10,
            ls="",
            marker="^",
            c=f"C{i}",
            ms=15,
            markeredgecolor="k",
            markeredgewidth=2,
        )
        ax.legend(prop={"size": 18}, loc="best")

    return ax
