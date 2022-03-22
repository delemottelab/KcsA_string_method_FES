"""This module contains functions to generate cvs."""

import numpy as np


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
    from numpy.linalg import norm

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
    from tqdm import tqdm

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
        d.append([ts.dt * i, norm(csel1 - csel2)/10])
    d = array(d)
    if not time:
        d = d[:, 1]
    return array(d)
