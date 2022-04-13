"""This module contains functions to generate cvs."""

import gc
import multiprocessing as mp
from glob import glob

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import dist
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


# class distance_pairs_av(AnalysisBase):
#     from numpy.linalg import norm

#     def __init__(self, u, mda_sel_txt0, mda_sel_txt1, tpr_path, verbose=False):
#         def wrap_and_prepare(u, trp_path):
#             import MDAnalysis.transformations as trans
#             from MDAnalysis.transformations.base import TransformationBase

#             # for reuniting all the subunits
#             class GroupHug(TransformationBase):
#                 def __init__(self, center, *others):
#                     super().__init__(max_threads=1, parallelizable=True)
#                     self.c = center
#                     self.o = others

#                 @staticmethod
#                 def calc_restoring_vec(ag1, ag2):
#                     box = ag1.dimensions[:3]
#                     dist = ag1.center_of_mass() - ag2.center_of_mass()

#                     return box * np.rint(dist / box)

#                 def _transform(self, ts):
#                     # loop over other atomgroups and shunt them into nearest image to center
#                     for i in self.o:
#                         rvec = self.calc_restoring_vec(self.c, i)
#                         i.translate(+rvec)
#                     return ts

#             u_bond = mda.Universe(tpr_path)
#             u.add_bonds(u_bond.bonds.to_indices())

#             u_prot = u.select_atoms("protein")
#             prot_chain_list = []
#             for chain in u_prot.segments:
#                 prot_chain_list.append(chain.atoms)
#             non_prot = u.select_atoms("not protein")

#             unwrap_trans = trans.unwrap(u.atoms)

#             prot_group = GroupHug(*prot_chain_list)
#             center_in_box = trans.center_in_box(u_prot)
#             wrap = trans.wrap(non_prot)
#             #        rot_fit_trans = trans.fit_rot_trans(u.select_atoms('name CA'), u.select_atoms('name CA'))
#             u.trajectory.add_transformations(
#                 *[unwrap_trans, center_in_box, wrap]  # prot_group,
#             )
#             return u

#         self.u = wrap_and_prepare(u.copy(), tpr_path)
#         super().__init__(trajectory=self.u.trajectory, verbose=verbose)
#         self.sel0 = self.u.select_atoms(mda_sel_txt0)
#         self.sel1 = self.u.select_atoms(mda_sel_txt1)
#         assert (
#             self.sel0.n_atoms == self.sel1.n_atoms
#         ), "Both atom selections should have the same number of atoms."
#         self.n_atoms = self.sel0.n_atoms
#         self.n_frames = u.trajectory.n_frames

#     def _prepare(self):
#         self.results.distances = np.zeros((self.n_frames, self.n_atoms))

#     def _single_frame(self):
#         def _dis(at0, at1):
#             return norm(at0.position - at1.position)

#         distances = []
#         for at0, at1 in zip(self.sel1.atoms, self.sel0.atoms):
#             distances.append(_dis(at0, at1))
#         self.results.distances[self._frame_index, :] = np.array(distances)

#     def _conclude(self):
#         self.results_pp = np.mean(self.results["distances"], axis=1)


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


class distance_pairs_av(AnalysisBase):
    def __init__(self, u, mda_sel_txt0, mda_sel_txt1, verbose=False):

        self.u = u
        super().__init__(trajectory=self.u.trajectory, verbose=verbose)
        self.sel0 = self.u.select_atoms(mda_sel_txt0)
        self.sel1 = self.u.select_atoms(mda_sel_txt1)
        assert (
            self.sel0.n_atoms == self.sel1.n_atoms
        ), "Both atom selections should have the same number of atoms."
        self.n_atoms = self.sel0.n_atoms
        self.n_frames = u.trajectory.n_frames

    def _prepare(self):

        self.results.distances = np.zeros((self.n_frames, self.n_atoms))

    def _single_frame(self):

        distances = dist(self.sel0, self.sel1, box=self._ts.dimensions)[2, :]
        self.results.distances[self._frame_index, :] = distances

    def _conclude(self):
        self.results_pp = np.mean(self.results["distances"], axis=1)
