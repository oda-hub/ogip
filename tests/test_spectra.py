import logging
import os
import numpy as np
import pytest

import ogip.spec
import ogip.core
import ogip.tools
from matplotlib import pylab as plt

logging.basicConfig(level=logging.DEBUG)


def test_read_something():

    for fn, c in [
        ("tests/data/phaI.fits.gz", ogip.spec.PHAI),
        ("tests/data/rmf_rt16_116.fits", ogip.spec.RMF),
        ("tests/data/rmf.fits.gz", ogip.spec.RMF),
        ("tests/data/arf.fits.gz", ogip.spec.ARF),
    ]:
        o = ogip.core.open_something(fn)
        assert isinstance(o, c)
        print(o.to_long_string())


def synth_e():
    return (np.logspace(1, 3, 100),
            np.logspace(1, 3, 50))


def synth_phaI():  # noqa: N802
    elh, eb = synth_e()

    r = np.random.normal(size=eb[:-1].shape) * 100

    phaI = ogip.spec.PHAI.from_arrays(  # noqa: N806
        1.,
        rate=r,
        stat_err=r / 100.
    )

    return phaI


def synth_rmf():
    elh, eb = synth_e()

    m_e, m_eb = np.meshgrid(elh[:-1], eb[:-1])

    d = (m_e - m_eb) / (0.1 * m_e)
    matrix = np.zeros_like(d)

    m = d < 10
    matrix[m] = np.exp(-d**2)[m]

    rmf = ogip.spec.RMF.from_arrays(
        energ_lo=elh[:-1],
        energ_hi=elh[1:],
        matrix=np.transpose(matrix),
        e_min=eb[:-1],
        e_max=eb[1:],
    )

    return rmf


def test_bad_arrays():
    pass


def test_from_arrays():

    rmf = synth_rmf()

    rmf.to_fits("rmf.fits")


@pytest.mark.skipif(not os.environ.get('HEADAS', None), reason="need HEASoft to test reading with xspec")
def test_read_xspec():
    import xspec

    synth_phaI().to_fits("phaI.fits")
    synth_rmf().to_fits("rmf.fits")

    s = xspec.Spectrum("phaI.fits")

    s.response = "rmf.fits"


def test_rebin():
    pha = ogip.core.open_something("tests/data/phaI.fits.gz")
    rmf = ogip.core.open_something("tests/data/rmf_rt16_116.fits")

    assert isinstance(pha, ogip.spec.PHAI)
    assert isinstance(rmf, ogip.spec.RMF)

    new_bins = ogip.spec.log_bins(5, 25, 250)
    rebinned_pha, rebinned_rmf = ogip.spec.rebin(pha, rmf, new_bins)

    assert rebinned_rmf._matrix.shape == (2466, 5)
    assert rebinned_pha._rate.shape == (5,)

    f = plt.figure()
    ogip.tools.plot(pha, lambda x: x, rmf, fig=f)
    ogip.tools.plot(rebinned_pha, lambda x: x, rebinned_rmf, fig=f)
    plt.savefig("png.png")


def test_unfolding():
    pha = ogip.core.open_something("tests/data/phaI.fits.gz")
    rmf = ogip.core.open_something("tests/data/rmf_rt16_116.fits")

    assert isinstance(pha, ogip.spec.PHAI)
    assert isinstance(rmf, ogip.spec.RMF)

    def model_gen(p):
        print("test", p)
        return lambda energy: (p[0] * (energy / 50)**p[1])

    ref_model = model_gen((1e-2, -2))

    pha._stat_err = pha._stat_err * 50
    pha._rate = ogip.tools.synthesise(pha, ref_model, rmf)

    mask = (rmf._e_min > 50) & (rmf._e_min < 100)

    r, f_model = ogip.tools.fit(model_gen, [0.5e-2, -3], [(pha, rmf, None, mask)])

    print(r, f_model)

    f = plt.figure()

    plt.plot(rmf._e_min, pha._rate / rmf.d_e)

    model_spec = ogip.tools.convolve(ref_model, rmf)
    plt.plot(rmf._e_min, model_spec / rmf.d_e, label="ref")

    model_spec = ogip.tools.convolve(f_model, rmf)
    plt.plot(rmf._e_min, model_spec / rmf.d_e, label="fitted")

    plt.plot(rmf._e_min, pha._rate / rmf.d_e)

    model_spec = ogip.tools.convolve(ref_model, rmf)
    plt.plot(rmf._e_min, model_spec / rmf.d_e)

    ll = ogip.tools.get_mloglike(pha, ref_model, rmf, mask=mask)

    plt.title(f"-loglike = {ll} {r}")

    plt.loglog()

    plt.savefig("png.png")

    plt.figure()

    ref_model(rmf._energ_lo)

    plt.plot(rmf._energ_lo, ref_model(rmf._energ_lo))

    plt.plot(rmf._e_min,
             pha._rate / model_spec * ref_model(rmf._e_min))

    plt.plot(rmf._e_min,
             pha._rate / rmf.d_e * ogip.tools.get_unfolding_factor(ref_model, rmf))

    plt.loglog()

    ogip.tools.plot(pha, ref_model, rmf, fig=f, unfolded=True)

    # TODO: check that it all looks the same

    plt.savefig("unf_png.png")


def test_read_poisson():
    import ogip.spec
    pha = ogip.spec.PHAI.from_file_name("tests/data/MOS1source_spectrum_150_rbn.pi")
    from astropy.io import fits as pf
    import numpy as np
    ff = pf.open("tests/data/MOS1source_spectrum_150_rbn.pi")
    counts = ff[1].data['COUNTS']

    assert(np.abs(np.sum(counts - (pha._rate)*pha._exposure)) < 1e-12)
    assert(np.sum(pha._rate*pha._exposure - (pha._stat_err*pha._exposure)**2))
