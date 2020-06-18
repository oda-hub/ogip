import logging
import os
import numpy as np
import pytest

import ogip.spec

logging.basicConfig(level=logging.DEBUG)

def test_read_something():
    import ogip.core

    for fn, c in [
                ("tests/data/phaI.fits.gz", ogip.spec.PHAI),
                ("tests/data/rmf_rt16_116.fits", ogip.spec.RMF),
                ("tests/data/rmf.fits.gz", ogip.spec.RMF),
            ]:
        o = ogip.core.open_something(fn)
        assert isinstance(o, c)
        print(o.to_long_string())



def synth_e():
    return (np.logspace(1, 3, 100),
            np.logspace(1, 3, 50))

def synth_phaI():
    elh, eb = synth_e()

    r = np.random.normal(size=eb[:-1].shape)*100

    phaI = ogip.spec.PHAI.from_arrays(
                1.,
                rate=r,
                stat_err=r/100.
            )

    return phaI

def synth_rmf():
    elh, eb = synth_e()

    m_e, m_eb = np.meshgrid(elh[:-1],eb[:-1])

    d = (m_e-m_eb)/(0.1*m_e)
    matrix = np.zeros_like(d)

    m = d<10
    matrix[m] = np.exp(-d**2)[m]

    rmf = ogip.spec.RMF.from_arrays(
                energ_lo = elh[:-1],
                energ_hi = elh[1:],
                matrix = np.transpose(matrix),
                e_min = eb[:-1],
                e_max = eb[1:],
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
