import logging
import ogip.lc
import os

logging.basicConfig(level=logging.DEBUG)


def test_read_write():
    files = ["tests/data/FPMA_3.0_7.0.flc",
             "tests/data/FPMA_7.0_30.0_sr.lc",
             "tests/data/IBIS_lc_Swift_J151857.0-572147_30-50.fits",
             "tests/data/JMX1_lc_Swift_J151857.0-572147_3-10.fits"
             ]

    for fn in files:
        lc = ogip.lc.Rate(fn)
        k = fn.rfind('.')
        new_fn = fn[0:k] + '_ogip.' + fn[k + 1:]
        lc.to_fits(new_fn)
        assert (os.path.isfile(new_fn))
