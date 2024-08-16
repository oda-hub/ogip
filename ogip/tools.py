# flake8: noqa

from typing import Callable, Optional, Tuple
import numpy as np

import emcee

import scipy
import scipy.stats
import scipy.optimize
from scipy.stats import norm

from ogip.spec import PHAI, RMF, ARF

import logging

logger = logging.getLogger(__name__)


def convolve(source_model: Callable, rmf: RMF, arf: Optional[ARF] = None, args=None):
    if arf is None:
        _arf = np.ones_like(rmf.c_energ)
    else:
        _arf = arf._arf

    if args is None:
        args = []

    r = (
        np.outer(
            _arf * source_model(rmf.c_energ, *args) * rmf.d_energ,
            np.ones_like(rmf._e_min),
        )
        * rmf._matrix
    ).sum(0)
    # r = (np.outer(_arf * source_model(rmf.c_energ)*rmf.d_energ, np.ones_like(rmf._e_min))*rmf._matrix).sum(0)

    logger.debug("convolved model (arf=%s): %s", arf, r)

    return r


def fit(
    source_model_generator: Callable,
    p0: list,
    spectra: list[tuple[PHAI, RMF, Optional[ARF], np.ndarray]],
    method="COBYLA",
):
    def F(p):
        r = 0

        for _ in spectra:
            r += get_mloglike(
                _[0], source_model_generator(p), rmf=_[1], arf=_[2], mask=_[3]
            )

        logger.debug("trying, result %s for parameters %s", r, p)
        # print(p)
        return r

    r = scipy.optimize.minimize(F, p0, method=method)

    logger.debug(
        "optimized model spec %s",
        convolve(source_model_generator(r.x), spectra[0][1], spectra[0][2]),
    )

    return r, source_model_generator(r.x)


def log_prior(p, p0):
    for requested_p, (pc, p1, p2) in zip(p, p0):
        if not (p1 < requested_p < p2):
            return -np.inf

    return 0.0


def log_prob(p, p0, spectra, source_model_with_pars):
    lp = log_prior(p, p0)
    if not np.isfinite(lp):
        return -np.inf
    else:
        r = 0
        for pha, rmf, arf, mask in spectra:
            convolved = convolve(source_model_with_pars, rmf, arf, args=(p,))

            if mask is None:
                mask = np.ones_like(convolved, dtype=bool)

            logprob = scipy.stats.norm(
                convolved, pha._stat_err).logpdf(pha._rate)

            for k in zip(convolved, pha._rate, pha._stat_err, logprob):
                logger.debug("convolved: %s rate %s err %s logprob %s", *k)

            r += np.sum(logprob[mask])

        if np.isnan(r):
            # print("this returns NaN:", p, r)
            return -np.inf
        else:
            # print("returning", r, lp)
            return r + lp


def sample(
    source_model_with_pars: Callable,
    p0: list[tuple[float, float, float]],
    spectra: list[tuple[PHAI, RMF, Optional[ARF], np.ndarray]],
    nwalkers=100,
    nsteps=10,
    n_processes=None,
):
    p0_c = np.array([[_[0] for _ in p0]] * nwalkers)
    p0_e = np.array([[(_[2] - _[1]) / 5 for _ in p0]] * nwalkers)
    initial_state = norm(p0_c, p0_e).rvs()
    print(initial_state.shape, len(p0))

    from multiprocessing import Pool

    def S(pool=None):
        sampler = emcee.EnsembleSampler(
            nwalkers,
            len(p0),
            log_prob,
            pool=pool,
            args=(p0, spectra, source_model_with_pars),
        )
        sampler.run_mcmc(initial_state, nsteps, progress=True)
        return sampler

    if n_processes == 1:
        sampler = S()
    else:
        with Pool(n_processes) as pool:
            sampler = S(pool)

    return sampler


def plot_chain(sampler, labels=None):
    from matplotlib import pylab as plt

    fig, axes = plt.subplots(sampler.ndim, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(sampler.ndim):
        try:
            ax = axes[i]
        except:
            ax = axes
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        if labels is not None:
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")


# def get_unfolding_factor(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF]=None) -> Tuple[np.ndarray, np.ndarray]:
def get_unfolding_factor(
    model: Callable, rmf: RMF, arf: Optional[ARF] = None
) -> Tuple[np.ndarray, np.ndarray]:
    model_spec = convolve(model, rmf, arf)

    e1 = rmf._e_min
    e2 = rmf._e_max

    return model(e1) / model_spec * rmf.d_e


def get_loglike(
    pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF] = None, mask=None
):
    convolved = convolve(model, rmf, arf)

    if mask is None:
        mask = np.ones_like(convolved, dtype=bool)

    logprob = scipy.stats.norm(convolved, pha._stat_err).logpdf(pha._rate)

    for k in zip(convolved, pha._rate, pha._stat_err, logprob):
        logger.debug("convolved: %s rate %s err %s logprob %s", *k)

    return np.sum(logprob[mask])


def get_mloglike(*args, **kwargs):
    return -1 * get_loglike(*args, **kwargs)


def synthesise(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF] = None):
    return scipy.stats.norm(convolve(model, rmf, arf), pha._stat_err).rvs()


def plot(
    pha: PHAI,
    model: Callable,
    rmf: RMF,
    arf: Optional[ARF] = None,
    fig=None,
    label_prefix=None,
    plot_kwargs={},
    unfolded=False,
    e_power=0,
    step=False,
    scale_factor=1,
    plot_model=False,
):
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pylab as plt

    ie1 = rmf._energ_lo
    ie2 = rmf._energ_hi
    e1 = rmf._e_min
    e2 = rmf._e_max

    model_spec = convolve(model, rmf, arf)
    logger.debug("model_spec:", model_spec)

    def prefix(l):
        if label_prefix is None:
            return l
        else:
            return f"{label_prefix}: {l}"

    if "label" not in plot_kwargs:
        plot_kwargs["label"] = prefix("model")

    if fig is None:
        plt.figure(figsize=(15, 10))

    if plot_model:
        if unfolded:
            plt.plot(
                rmf.c_energ, model(rmf.c_energ) *
                rmf.c_energ**e_power * scale_factor
            )
        else:
            # print("model_spec, rmf.d_e, scale_factor", model_spec, rmf.d_e, scale_factor)
            # print("model_spec/rmf.d_e * scale_factor", model_spec/rmf.d_e * scale_factor)

            plt.step(
                e1, model_spec / rmf.d_e * scale_factor, where="post", **plot_kwargs
            )

    if pha is not None:
        if unfolded:
            unfolding_factor = get_unfolding_factor(model, rmf, arf)
        else:
            unfolding_factor = 1

        if step:
            x = plt.step(
                rmf.c_e,
                pha._rate
                / rmf.d_e
                * unfolding_factor
                * rmf.c_e**e_power
                * scale_factor,
                where="mid",
                **plot_kwargs,
            )
            color = x[0].get_color()
        else:
            color = None

        rate_or_uplim = (
            pha._rate / rmf.d_e * unfolding_factor * rmf.c_e**e_power * scale_factor
        )
        err = (
            pha._stat_err
            / rmf.d_e
            * unfolding_factor
            * rmf.c_e**e_power
            * scale_factor
        )

        uplims = rate_or_uplim - err < 0
        rate_or_uplim[uplims] = err[uplims] * 3  # factor!

        plt.errorbar(
            rmf.c_e,
            rate_or_uplim,
            err,
            xerr=rmf.d_e / 2,
            uplims=uplims,
            ls="",
            c=color,
            **plot_kwargs,
        )

    plt.legend()
    plt.loglog()

    plt.xlim([15, 600])
    plt.grid()


def transform_rmf(
    rmf: RMF,
    arf: ARF,
    bias_function: Callable,
    preserved_projection=None,
    tolerance=0.05,
) -> RMF:
    new_rmf = RMF.from_arrays(
        energ_lo=bias_function(rmf._energ_lo.copy()),
        energ_hi=bias_function(rmf._energ_hi.copy()),
        e_min=rmf._e_min.copy(),
        e_max=rmf._e_max.copy(),
        matrix=rmf._matrix.copy(),
    )

    if preserved_projection is not None:
        # projection on model should be preserved
        corr = convolve(preserved_projection, new_rmf, arf) / convolve(
            preserved_projection, rmf, arf
        )

        new_rmf._matrix /= np.outer(np.ones_like(new_rmf._energ_lo), corr)

        s = convolve(preserved_projection, rmf, arf)
        t_s = convolve(preserved_projection, new_rmf, arf)

        assert np.all(np.abs((s - t_s) / s) < tolerance)

    return new_rmf


# Roques & Jourdain 2018
def crab_ph_cm2_s_kev(en):
    K = 7.417e-4
    al1 = -1.98
    al2 = -2.33
    Ec = 500.0
    f = K * (en / 100) ** al1 * (np.exp(-en / Ec))
    m = en > Ec * (al1 - al2)
    f[m] = (
        K
        * ((al1 - al2) * Ec / 100) ** (al1 - al2)
        * (en / 100) ** al2
        * np.exp(-(al1 - al2))
    )[m]

    return f


def le_bias_function(energy, bias_keV):
    # TODO: preserves 60 and 511, shifts LE
    return energy + bias_keV
