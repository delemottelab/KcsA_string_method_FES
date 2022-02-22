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


def path_cvs(coord, path, lam=None):

    return path_cvs
