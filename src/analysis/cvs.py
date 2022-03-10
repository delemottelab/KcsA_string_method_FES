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
