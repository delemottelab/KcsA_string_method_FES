import multiprocessing as mp
import sys

import numpy as np
from stringmethod.config import *
from stringmethod.postprocessing import *


def load_swarm_data(extract, first_iteration=1, last_iteration=None):
    if last_iteration is None:
        last_iteration = sys.maxsize
    if extract:
        config = load_config("config.json")

        ce = CvValueExtractor.from_config(
            config=config,
            # Exclude the first iterations to let the system equilibrate.
            first_iteration=first_iteration,
            # Usefull to make blocks of the simulation
            last_iteration=last_iteration,
        )
        ce.run()
        ce.persist()
    return np.load("postprocessing/cv_coordinates.npy")


def k_means_cluster(data, k, stride=1, max_iter=500, n_jobs=1, seed=None):
    from deeptime.clustering import KMeans

    estimator = KMeans(
        n_clusters=k,  # place 100 cluster centers
        init_strategy="uniform",  # uniform initialization strategy
        max_iter=max_iter,  # don't actually perform the optimization, just place centers
        fixed_seed=seed,
        n_jobs=n_jobs,
    )
    clusters = (
        estimator.fit(data[::stride, :, :].reshape(-1, data.shape[2]))
        .fetch_model()
        .transform(data.reshape(-1, data.shape[2]))
    )
    clusters = clusters.reshape(-1, data.shape[1])
    return clusters


def get_vamp_vs_k(
    n_clustercenters,
    data,
    n_jobs,
    allow_failed_msms=False,
    reversible=True,
    scores=None,
):
    import logging

    import deeptime.markov as markov
    import matplotlib.pyplot as plt
    from deeptime.decomposition import vamp_score_cv
    from deeptime.util import confidence_interval
    from tqdm.autonotebook import tqdm

    plt.rcParams["axes.facecolor"] = "#f9f9fb"
    plt.rcParams["grid.color"] = "white"
    plt.rcParams["grid.linestyle"] = "-"
    plt.rcParams["grid.linewidth"] = 2
    plt.rcParams["axes.grid"] = True
    plt.rcParams["lines.solid_capstyle"] = "round"

    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        logger.setLevel(logging.ERROR)

    n_iter = 5
    if scores is None:
        scores = np.zeros((len(n_clustercenters), n_iter))
        for n, k in tqdm(
            enumerate(n_clustercenters),
            total=len(n_clustercenters),
            desc="Loop over k:",
        ):
            failed_counter = 0
            for m in tqdm(range(n_iter), desc="Loop over iterations", leave=False):
                _cl = k_means_cluster(data, k, stride=10, max_iter=50, n_jobs=n_jobs)

                # TODO, refactor
                if allow_failed_msms:
                    try:
                        estimator = markov.msm.MaximumLikelihoodMSM(
                            reversible=reversible,
                            stationary_distribution_constraint=None,
                            lagtime=1,
                        )

                        counts = (
                            markov.TransitionCountEstimator(
                                lagtime=1, count_mode="sample"
                            )
                            .fit(_cl, n_jobs=n_jobs)
                            .fetch_model()
                        )
                        _msm = estimator.fit(counts, n_jobs=n_jobs)
                        # return _msm, _cl
                        # exit

                        scores[n, m] = vamp_score_cv(
                            _msm,
                            trajs=[c for c in _cl],
                            n=1,
                            lagtime=1,
                            dim=min(10, k),
                            n_jobs=n_jobs,
                        )[0]
                    except:
                        failed_counter += 1
                        scores[n, m] = np.nan
                else:
                    estimator = markov.msm.MaximumLikelihoodMSM(
                        reversible=reversible,
                        stationary_distribution_constraint=None,
                        lagtime=1,
                    )

                    counts = (
                        markov.TransitionCountEstimator(lagtime=1, count_mode="sample")
                        .fit(_cl, n_jobs=n_jobs)
                        .fetch_model()
                    )
                    _msm = estimator.fit(counts, n_jobs=n_jobs)
                    # return _msm, _cl
                    # exit

                    scores[n, m] = vamp_score_cv(
                        _msm,
                        trajs=[c for c in _cl],
                        n=1,
                        lagtime=1,
                        dim=min(10, k),
                        n_jobs=n_jobs,
                    )[0]
            if allow_failed_msms:
                print(f"For k={k}, {failed_counter} iterations failed out of {n_iter}.")

        # Plotting
    fig, ax = plt.subplots(1, 1)
    lower, upper = confidence_interval(scores.T.tolist(), conf=0.9)
    ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
    ax.plot(n_clustercenters, np.nanmean(scores, axis=1), "-o")
    ax.semilogx()
    ax.set_xlabel("number of cluster centers")
    ax.set_ylabel("VAMP-2 score")
    fig.tight_layout()

    return fig, ax, scores


def get_bayesian_msm(clusters, n_samples=100):
    import deeptime.markov as markov

    estimator = markov.msm.BayesianMSM(
        n_samples=100, reversible=True, stationary_distribution_constraint=None
    )

    counts = (
        markov.TransitionCountEstimator(lagtime=1, count_mode="effective")
        .fit(clusters)
        .fetch_model()
    )

    msm = estimator.fit(counts).fetch_model()
    weights = msm.gather_stats(
        "compute_trajectory_weights", dtrajs=clusters.reshape(-1, 1)
    )

    return msm, weights


def get_msm(clusters, n_jobs=1, reversible=True):
    import logging

    import deeptime.markov as markov

    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        logger.setLevel(logging.ERROR)

    estimator = markov.msm.MaximumLikelihoodMSM(
        reversible=reversible,
        # stationary_distribution_constraint=None,
        lagtime=1,
    )
    counts = markov.TransitionCountEstimator.count(
        lagtime=1, count_mode="sample", dtrajs=clusters, n_jobs=n_jobs
    )
    counts = counts.astype(int)
    msm = estimator.fit(counts, n_jobs=n_jobs).fetch_model()
    weights = np.array(
        msm.compute_trajectory_weights(np.concatenate(clusters.reshape(-1, 1)))
    ).squeeze()

    return msm, weights


def get_kde(samples, weights=None, bandwidth=None, nbins=55, extent=None):
    assert samples.shape[-1] < 3, "Can only do 1D or 2D KDE"
    if samples.shape[-1] == 1:
        return get_kde_1d(
            samples, weights=weights, bandwidth=bandwidth, nbins=nbins, extent=extent
        )
    elif samples.shape[-1] == 2:
        return get_kde_2d(
            samples, weights=weights, bandwidth=bandwidth, nbins=nbins, extent=extent
        )


def get_kde_2d(samples, weights=None, bandwidth=None, nbins=55, extent=None):
    """
    aklsjfa;sldjfas
    """

    from scipy.stats import gaussian_kde

    samples = samples.reshape(-1, samples.shape[1])
    if extent is None:
        xmin = samples[:, 0].min()
        xmax = samples[:, 0].max()
        ymin = samples[:, 1].min()
        ymax = samples[:, 1].max()
    else:
        xmin, xmax, ymin, ymax = extent
    nbins = nbins * 1j
    X, Y = np.mgrid[xmin:xmax:nbins, ymin:ymax:nbins]
    positions = np.vstack([X.ravel(), Y.ravel()])
    kernel = gaussian_kde(
        samples.T,
        weights=weights,
        bw_method=bandwidth,
    )
    Z = np.reshape(kernel(positions), X.shape)
    Z = Z.T
    extent = [xmin, xmax, ymin, ymax]

    return Z, extent


def get_kde_2d_custom(
    samples, weights=None, bandwidth=None, nbins=55, extent=None, density=True, progressbar=True
):
    """
    aklsjfa;sldjfas
    """

    from scipy.stats import multivariate_normal
    from tqdm.autonotebook import tqdm

    samples = samples.reshape(-1, 2)
    if extent is None:
        xmin = samples[:, 0].min()
        xmax = samples[:, 0].max()
        ymin = samples[:, 1].min()
        ymax = samples[:, 1].max()
    else:
        xmin, xmax, ymin, ymax = extent
    X, Y = np.mgrid[xmin : xmax : nbins * 1j, ymin : ymax : nbins * 1j]
    positions = np.dstack((X, Y))
    Z = np.zeros([nbins, nbins])

    assert np.isscalar(bandwidth) or bandwidth.shape == (
        2,
        2,
    ), "The bandwidth must be a scalr or a 2x2 covariance matrix."
    if np.isscalar(bandwidth):
        cov = np.cov(samples.T, aweights=weights, bias=False)
        cov = cov * bandwidth**2
    else:
        cov = bandwidth

    for w, sample in tqdm(
        zip(weights, samples), total=weights.shape[0], desc="Loop over samples", disable=not progressbar
    ):
        rv = multivariate_normal([sample[0], sample[1]], cov)
        Z += w * rv.pdf(positions)

    extent = [xmin, xmax, ymin, ymax]
    if density:
        norm = np.sum(Z)
        dx = (xmax - xmin) / (nbins - 1)
        dy = (ymax - ymin) / (nbins - 1)
        norm = norm * dx * dy
        Z = Z / norm
    return Z.T, extent, cov


def get_kde_1d(samples, weights=None, bandwidth=None, nbins=55, extent=None):
    """
    aklsjfa;sldjfas
    """

    from scipy.stats import gaussian_kde

    samples = samples.flatten()
    if extent is None:
        xmin = samples.min()
        xmax = samples.max()
    else:
        xmin, xmax = extent
    positions = np.linspace(xmin, xmax, nbins)
    kernel = gaussian_kde(
        samples,
        weights=weights,
        bw_method=bandwidth,
    )
    extent = [xmin, xmax]
    Z = kernel(positions)

    return Z, extent


def _get_bootstrap(n_boot, block_length, clusters, cv_proj, bandwidth, extent, nbin):

    random = np.random.choice(n_boot, n_boot)
    mask = []
    for r in random:
        mask += list(range(r * block_length, (r + 1) * block_length))
    try:
        # print("do msm")
        _, w = get_msm(clusters[mask])
    except:
        # print("failed msm")
        h = np.zeros([nbin, nbin])
        h = np.nan
        error = 1
        return h, error
    # print("do KDE ")
    h, _, = get_kde(
        cv_proj[mask, :, :],
        w,
        bandwidth,
        extent=extent,
        nbins=nbin,
    )
    error = 0

    return h, error


def get_error(
    cv_proj,
    clusters,
    extent,
    n_boot=200,
    bandwidth=0.05,
    nbin=55,
    blocks=[2, 4, 8, 16, 32, 64],
    n_jobs=1,
    seed=None,
):
    from functools import partial

    from tqdm.autonotebook import tqdm

    if seed is not None:
        np.random.seed(seed)

    assert (
        type(bandwidth) == float
    ), "You have to choose a bandwidth, otherwise the bandwith will change between bootstrap iterations"
    ndat = cv_proj.shape[0]
    errors = []
    for b in tqdm(blocks, desc="Loop over blocks"):
        block_length = ndat // b
        get_bootstrap = partial(
            _get_bootstrap, b, block_length, clusters, cv_proj, bandwidth, extent, nbin
        )

        hist_err = []
        with mp.Pool(n_jobs) as tp:
            hist_err = tp.starmap(
                get_bootstrap,
                tqdm([() for _ in np.arange(n_boot)], total=n_boot, desc="bootstraps"),
            )
        # print("done looping")
        # hist_err = []
        # for _ in tqdm(np.arange(n_boot), total=n_boot, desc="bootstraps"):
        #     hist_err.append(get_bootstrap())
        err = 0
        histograms = []
        for h, e in hist_err:
            err += e
            histograms.append(h)
        print(f"{err/n_boot *100} of MSMs failed")
        # return histograms
        histograms = np.array(histograms)
        x_mean = np.nanmean(histograms, axis=0)
        std_err = np.nanstd(histograms, axis=0, ddof=1)
        # Dividing by x_mean propagates the uncertainty from histogram uncertainty to
        # free energy uncertainty.
        std_err = std_err / x_mean

        errors.append(std_err)
    errors = np.array(errors)
    errors[errors == np.inf] = np.nan

    return errors


def get_hdi(x, axis, alpha=0.06):
    x_mean = np.nanmedian(x, axis=axis)
    percentiles = 100 * np.array([alpha / 2.0, 1.0 - alpha / 2.0])
    hdi = np.nanpercentile(x, percentiles, axis=axis)
    return x_mean, hdi


def project_property_on_cv_kde(
    cv_proj, weights, proper, bandwidth=None, nbins=55, normalize=True, F_cutoff_KT=40, progressbar=True
):

    n_data = cv_proj.shape[0] * cv_proj.shape[1]
    if len(proper.shape) == 3:
        proper = proper.reshape(n_data)
    elif len(proper.shape) != 1:
        raise Exception(
            "The shape of property should be (ndata,) or (n_iter, n_frames_per_iter, 1)"
        )

    weights_proper = proper * weights

    count_of_cv, extent, cov = get_kde_2d_custom(
        cv_proj, weights, bandwidth, nbins=nbins, density=False, progressbar=progressbar
    )
    prop_of_cv, extent, cov = get_kde_2d_custom(
        cv_proj, weights_proper, cov, nbins=nbins, density=False, progressbar=progressbar
    )
    if normalize:
        prop_of_cv[count_of_cv < np.exp(-F_cutoff_KT)] = np.nan
        prop_of_cv = prop_of_cv / count_of_cv

    return prop_of_cv, extent
