import logging
import numpy as np # type: ignore
import astropy.io.fits as fits # type: ignore

logger = logging.getLogger()

from ogip.core import Reader

class Spectrum(Reader):

    @staticmethod
    def from_file_name(fn):
        c = Spectrum()
        return c

    def to_long_string(self):
        return repr(self)

    ###

class PHAI(Spectrum):
    @staticmethod
    def from_arrays(exposure, rate=None, stat_err=None, sys_err=None, quality=None, counts=None):
        self = PHAI()

        self._exposure = exposure

        if counts is not None:
            logger.warning("counts found: converting to rate")
            rate = counts/exposure

            if stat_err is None:
                _stat_err = counts**0.5
                logger.warning("assuming poisson stat errors. Is it really what you want?")

        if rate is None:
            raise Exception("need rate or counts")

        self._rate = rate

        if stat_err is not None:
            self._stat_err = stat_err
        else:
            raise Exception("need stat err")

        if sys_err is not None:
            self._sys_err = sys_err
        else:
            logger.warning("assuming no systematic errors. Probably ok")
            self._sys_err = np.zeros_like(rate)

        if quality is not None:
            self._quality = quality
        else:
            self._quality = np.ones_like(rate)

        return self

    @property
    def spectrum_hdu(self):

        return fits.BinTableHDU.from_columns([
                    fits.Column(name='RATE', array=self._rate, format='1E'),
                    fits.Column(name='STAT_ERR', array=self._stat_err, format='1E'),
                    fits.Column(name='SYS_ERR', array=self._sys_err, format='1E'),
                    fits.Column(name='QUALITY', array=self._quality, format='1I'),
                ],
                header=fits.Header(cards=dict(
                                      EXTNAME="SPECTRUM", # - the name (i.e. type) of the extension
                                      TELESCOP="", # - the "telescope" (i.e. mission/satellite name).
                                      INSTRUME="", #- the instrument/detector.
                                      #FILTER - the instrument filter in use (if any)
                                      AREASCAL=1.,
                                      BACKSCAL=1.,
                                      EXPOSURE=self._exposure, # the integration time (in seconds) for the PHA data (assumed to be corrected for deadtime, data drop-outs etc. )
                                      BACKFILE="", #- the name of the corresponding background file (if any)
                                      #CORRFILE - the name of the corresponding correction file (if any)
                                      #CORRSCAL - the correction scaling factor.
                                      RESPFILE="", # - the name of the corresponding (default) redistribution matrix file (RMF; see George et al. 1992a).
                                      ANCRFILE="", # - the name of the corresponding (default) ancillary response file (ARF; see George et al. 1992a).
                                      HDUCLASS="OGIP", # - should contain the string "OGIP" to indicate that this is an OGIP style file.
                                      HDUCLAS1="SPECTRUM", # - should contain the string "SPECTRUM" to indicate this is a spectrum.
                                      HDUVERS="1.2.1", #- the version number of the format (this document describes version 1.2.1)
                                      #POISSERR #- whether Poissonian errors are appropriate to the data (see below).
                                      #CHANTYPE #- whether the channels used in the file have been corrected in anyway (see below).
                                      DETCHANS=len(self._rate), #- the total number of detector channels available.
                                  )),
            )

    def to_fits(self, fn: str):
        logger.info("store to fits %s %s", self, fn)

        fits.HDUList([
                fits.PrimaryHDU(),
                self.spectrum_hdu,
            ]).writeto(fn, overwrite=True)


class PHAII(Spectrum):
    pass


class RMF:

    def __init__(self):
        pass

    @staticmethod
    def from_arrays(energ_lo: float, energ_hi, matrix, e_min, e_max):
        self = RMF()

        self._energ_lo = energ_lo # type: ignore
        self._energ_hi = energ_hi # type: ignore
        self._matrix = matrix # type: ignore
        self._e_min = e_min # type: ignore
        self._e_max = e_max # type: ignore

        return self

    @property
    def matrix_hdu(self):
        return fits.BinTableHDU.from_columns([
                    fits.Column(name='ENERG_LO', array=self._energ_lo, format='1E'),
                    fits.Column(name='ENERG_HI', array=self._energ_hi, format='1E'),
                ],
                header=fits.Header(cards=dict(
                                      EXTNAME="MATRIX"
                                  )),
            )
    
    @property
    def ebounds_hdu(self):
        return fits.BinTableHDU.from_columns([
                    fits.Column(name='E_MIN', array=self._e_min, format='1E'),
                    fits.Column(name='E_MAX', array=self._e_max, format='1E'),
                ],
                header=fits.Header(cards=dict(
                                      EXTNAME="MATRIX"
                                  )),
            )
    
    def to_fits(self, fn):
        logger.info("store to fits %s %s", self, fn)

        fits.HDUList([
                fits.PrimaryHDU(),
                self.ebounds_hdu,
                self.matrix_hdu,
            ]).writeto(fn, overwrite=True)

class ARF:
    pass
