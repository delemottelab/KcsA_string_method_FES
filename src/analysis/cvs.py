"""This module contains functions to generate cvs."""

import gc
import multiprocessing as mp
from glob import glob

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from numpy.linalg import norm
from tqdm.autonotebook import tqdm

from .utils import natural_sort


def strings_to_SF_IG(strings, SF_cv_numbers, IG_cv_numbers):

    strings_to_SF_IG = np.hstack(
        [
            np.mean(strings[:, SF_cv_numbers, :], axis=1, keepdims=True),
            np.mean(strings[:, IG_cv_numbers, :], axis=1, keepdims=True),
        ]
    )
    return strings_to_SF_IG


def cvs_to_SF_IG(cv_coordinates, SF_cv_numbers, IG_cv_numbers):

    cvs_to_SF_IG = np.concatenate(
        [
            np.mean(
                [
                    cv_coordinates[:, :, SF_cv_numbers[0] : SF_cv_numbers[0] + 1],
                    cv_coordinates[:, :, SF_cv_numbers[1] : SF_cv_numbers[1] + 1],
                ],
                axis=0,
            ),
            np.mean(
                [
                    cv_coordinates[:, :, IG_cv_numbers[0] : IG_cv_numbers[0] + 1],
                    cv_coordinates[:, :, IG_cv_numbers[1] : IG_cv_numbers[1] + 1],
                ],
                axis=0,
            ),
        ],
        axis=2,
    )

    return cvs_to_SF_IG


def MSD_metric(v1, v2):

    return np.sum((v1 - v2) ** 2)


def get_path_lambda(path, metric=MSD_metric):

    n_beads = path.shape[1]
    lam = (
        2.3
        * (n_beads - 1)
        / np.sum([metric(path[:, i], path[:, i + 1]) for i in range(0, n_beads - 1)])
    )
    return lam


def cvs_to_path(vec, path, lam, metric=MSD_metric):
    from numpy.linalg import norm

    n_beads = path.shape[1]

    # The formula of the Branduardi paper indeces from i=1 to P, if j=i-1, then j=0 to P-1
    array = np.array([np.exp(-lam * metric(vec, bead)) for bead in path[:, :].T])

    s = np.sum(np.arange(1, n_beads + 1) * array) / np.sum(array) / (n_beads - 1)
    z = -np.log(np.sum(array)) / lam
    return np.array([s, z])


def average_strings_to_cv(cv_coordinates, cv_id_list):
    cv_id_list0 = []
    if type(cv_id_list[0]) is not list:
        cv_id_list0.append(cv_id_list)
    av = np.concatenate(
        [
            np.mean([cv_coordinates[:, :, c : c + 1] for c in cv], axis=0)
            for cv in cv_id_list0
        ],
        axis=2,
    )
    return av


def distance_atom_groups(
    u, sel1, sel2, time=False, progressbar=False, center_of_mass=False
):
    """
    Calculate the distance between the centers of geometry (or mass) between two groups (sel1, sel2) as a function of time in the trajectory trj.

    Parameters
    ----------
    u: MDA universe to analyz trajectory to analyze.
    sel1: MDA selection containing at least 1 atom.
    sel2: MDA selection containing at least 1 atom.
    center_of_mass: Use the center of mass instead of center of geometry.
    progressbar: Show progressbar.

    Returns
    -------
    d: numpy array.
    """
    from MDAnalysis import AtomGroup, Universe
    from numpy import array
    from numpy.linalg import norm

    assert isinstance(u, Universe), "u should be a MDAnlaysis universe."
    assert isinstance(sel1, AtomGroup), "sel1 should be a MDAnlaysis universe."
    assert isinstance(sel2, AtomGroup), "sel2 should be a MDAnlaysis universe."
    assert isinstance(progressbar, bool), "progressbar should be boolean."
    assert sel1.n_atoms >= 1, "sel1 should have at least 1 atom."
    assert sel2.n_atoms >= 1, "sel2 should have at least 1 atom."

    d = []
    for i, ts in tqdm(
        enumerate(u.trajectory), total=u.trajectory.n_frames, disable=not progressbar
    ):
        if center_of_mass:
            csel1 = sel1.center_of_mass()
            csel2 = sel2.center_of_mass()
        else:
            csel1 = sel1.centroid()
            csel2 = sel2.centroid()
        d.append([ts.dt * i, norm(csel1 - csel2) / 10])
    d = array(d)
    if not time:
        d = d[:, 1]
    return array(d)


class janin_chi1_av:
    def __init__(self, u, **kwargs):
        from MDAnalysis.analysis.dihedrals import Janin

        self.janin_res = u.select_atoms(kwargs["mda_sel_txt"])
        self.janin = Janin(self.janin_res, **kwargs)

    def run(self):
        self.janin.run()
        self.results_pp = np.mean(self.janin.results["angles"][:, :, 0], axis=1)

    def plot(self):
        assert (
            self.results_pp is not None
        ), "Need to run and postprocess before plotting."
        plt.plot(self.results_pp)


def _get_swarm(p, path_topology, mda_object, n_swarms, n_beads, **kwargs):
    swarms = []
    for b in range(1, n_beads - 1):
        for s in range(n_swarms):
            swarms.append(f"{p}/{b}/s{s}/traj_comp.xtc")
    # print(swarms)
    iter_num = p.split("/")[-2]
    try:
        u = mda.Universe(path_topology, swarms, verbose=True, refresh_offsets=True)
        if u.trajectory.n_frames != n_swarms * 2 * (n_beads - 2):
            print(f"Iteration {p} seems to be missing frames.")
    except IOError:
        print(f"Problem with iteration: {iter_num}")
        return None

    obj = mda_object(u, **kwargs)
    obj.run()
    results = obj.results_pp
    del u
    gc.collect()
    return results


def loop_over_iter(
    path,
    path_topology,
    mda_object,
    start=1,
    stop=None,
    step=1,
    n_jobs=1,
    **kwargs,
):

    from functools import partial

    paths = natural_sort(glob(f"{path}/md/*/"))
    if stop is None:
        stop = len(paths)
    paths = paths[start:stop:step]
    # print(paths)
    get_swarm = partial(
        _get_swarm,
        path_topology=path_topology,
        mda_object=mda_object,
        n_swarms=32,
        n_beads=18,
        **kwargs,
    )

    if n_jobs > 1:
        ## PARALLEL
        with mp.Pool(n_jobs) as tp:
            data = list(
                tqdm(
                    tp.imap(
                        get_swarm,
                        paths,
                    ),
                    total=len(paths),
                    desc="Swarms",
                )
            )
        ##
    else:
        ## SERIAL
        data = []
        for p in tqdm(paths, total=len(paths), desc="Swarms"):
            data.append(get_swarm(p))
        ##

    data = [d for d in data if d is not None]
    data = np.hstack(data)
    data = data.reshape([data.shape[0], 1])
    return data
