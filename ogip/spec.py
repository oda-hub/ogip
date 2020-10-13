import logging
import numpy as np # type: ignore
import astropy.io.fits as fits # type: ignore

logger = logging.getLogger()


class Spectrum:

    @staticmethod
    def from_file_name(fn):
        raise NotImplementedError

    def to_long_string(self):
        return repr(self)

    ###

class Readable:
    @staticmethod
    def from_file_name(cls, fn):
        for n, m in cls.__dict__.items():
            if n.startswith('from_file_name_'):
                try:
                    return m.__func__(fn)
                except Exception as e:
                    logger.debug("failed to read with %s.%s %s: %s", cls.__name__, n, m, e)

        raise Exception("failed to read with any method")

class PHAI(Spectrum):

    @staticmethod
    def from_file_name(fn):
        return Readable.from_file_name(PHAI, fn)

    @staticmethod
    def from_file_name_osa(fn):
        f = fits.open(fn)
        for e in f:
            try:
                return PHAI.from_arrays(
                            exposure=e.header['EXPOSURE'],
                            rate=e.data['RATE'],
                            stat_err=e.data['STAT_ERR'],
                            sys_err=e.data['SYS_ERR'],
                        )
            except Exception as ex:
                logger.debug("failed to read from %s: %s", e, ex)

        raise Exception("unable to read from any extension")

    @staticmethod
    def from_file_name_normal(fn):
        f = fits.open(fn)
        return PHAI.from_arrays(
                    exposure=f['SPECTRUM'].header['EXPOSURE'],
                    rate=f['SPECTRUM'].data['RATE'],
                    stat_err=f['SPECTRUM'].data['STAT_ERR'],
                    sys_err=f['SPECTRUM'].data['SYS_ERR'],
                )
    


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
            self._quality = np.zeros_like(rate)

        return self

    @property
    def spectrum_hdu(self):
        return fits.BinTableHDU.from_columns([
                    fits.Column(name='RATE', array=self._rate, format='1E'),
                    fits.Column(name='STAT_ERR', array=self._stat_err, format='1E'),
                    fits.Column(name='SYS_ERR', array=self._sys_err, format='1E'),
                    fits.Column(name='QUALITY', array=self._quality, format='1I'),
                ],
                # https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node6.html
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

    @property
    def total_counts(self):
        return np.nansum(self._rate)

    def to_long_string(self):
        return f"{self.__class__.__name__}: exposure {self._exposure} s, total count/s {self.total_counts}, {len(self._rate)} channels "


class PHAII(Spectrum):
    pass


class RMF:
    _telescop="not-a-telescope"
    _instrume="not-an-instrument"
    
    def to_long_string(self):
        return f"{self.__class__.__name__}: {len(self._e_min)} x {len(self._energ_lo)}  channels "

    @staticmethod
    def from_file_name(fn):
        return Readable.from_file_name(RMF, fn)
    
    @staticmethod
    def from_file_name_osaisgri(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
                    energ_lo=f['ISGR-RMF.-RSP'].data['ENERG_LO'],
                    energ_hi=f['ISGR-RMF.-RSP'].data['ENERG_HI'],
                    matrix=f['ISGR-RMF.-RSP'].data['MATRIX'],
                    e_min=f['ISGR-EBDS-MOD'].data['E_MIN'],
                    e_max=f['ISGR-EBDS-MOD'].data['E_MAX'],
                )
    
    @staticmethod
    def from_file_name_normal(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
                    energ_lo=f['MATRIX'].data['ENERG_LO'],
                    energ_hi=f['MATRIX'].data['ENERG_HI'],
                    matrix=np.vstack(f['MATRIX'].data['MATRIX']),
                    e_min=f['EBOUNDS'].data['E_MIN'],
                    e_max=f['EBOUNDS'].data['E_MAX'],
                )

    @staticmethod
    def from_arrays(energ_lo, energ_hi, matrix, e_min, e_max):
        self = RMF()


        if not (len(energ_lo) == len(energ_hi) == matrix.shape[0]):
            raise Exception(f"incompatible dimensions of mc energy, bounds {len(energ_lo)} {len(energ_hi)} but matrix {matrix.shape[0]}!")
        
        if not (len(e_min) == len(e_max) == matrix.shape[1]):
            raise Exception(f"incompatible dimensions of channel energy, bounds {len(e_min)} {len(e_max)} but matrix {matrix.shape[1]}!")

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
                    fits.Column(name='N_GRP', array=np.ones_like(self._energ_lo), format='1I'),
                    fits.Column(name='F_CHAN', array=0*np.ones_like(self._energ_lo), format='1I'),
                    fits.Column(name='N_CHAN', array=len(self._e_min)*np.ones_like(self._energ_lo), format='1I'),
                    fits.Column(name='MATRIX', array=self._matrix, format='PE'),
                    #fits.Column(name='MATRIX', array=self._matrix, format=f'{len(self._e_min)}E'),
                ],
                header=fits.Header(cards=dict(
                    # https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc7.1.1
                                      EXTNAME = 'MATRIX', #           / name of this binary table extension
                                      #BITPIX  =                    8 / 8-bit bytes
                                      #NAXIS   =                    2 / 2-dimensional binary table
                                      #NAXIS1  =                   34 / width of table in bytes
                                      #NAXIS2  =                 1180 / number of rows in table
                                      #PCOUNT  =              1031160 / Number of bytes acumulated in heap
                                      #GCOUNT  =                    1 / one data group (required keyword)
                                      #TFIELDS =                    6 / number of fields in each row
                                      #TTYPE1  = 'ENERG_LO'           / label for field   1
                                      #TFORM1  = 'E       '           / data format of field: 4-byte REAL
                                      #TUNIT1  = 'keV     '           / physical unit of field
                                      #TTYPE2  = 'ENERG_HI'           / label for field   2
                                      #TFORM2  = 'E       '           / data format of field: 4-byte REAL
                                      #TUNIT2  = 'keV     '           / physical unit of field
                                      #TTYPE3  = 'N_GRP   '           / label for field   3
                                      #TFORM3  = 'I       '           / data format of field: 2-byte INTEGER
                                      #TTYPE4  = 'F_CHAN  '           / label for field   4
                                      #TFORM4  = 'PI(2)   '           / data format of field: variable length array
                                      #TTYPE5  = 'N_CHAN  '           / label for field   5
                                      #TFORM5  = 'PI(2)   '           / data format of field: variable length array
                                      #TTYPE6  = 'MATRIX  '           / label for field   6
                                      #TFORM6  = 'PE(418) '           / data format of field: variable length array
                                      TLMIN4=0, #/ First legal channel number
                                      TLMAX4=len(self._e_min), #/ Highest legal channel number
                                      TELESCOP=self._telescop, #          / mission/satellite name
                                      INSTRUME=self._instrume, #           / instrument/detector
                                      #FILTER  = 'NONE    '           / filter information
                                      CHANTYPE='PI', #           / Type of channels (PHA, PI etc)
                                      DETCHANS=len(self._e_min), # / Total number of detector PHA channels
                                      LO_THRES=1.00E-07, # / Lower probability density threshold for matrix
                                      HDUCLASS='OGIP', #/KeywordinformationforCaltoolsSoftware.
                                      HDUCLAS1='RESPONSE', #/KeywordinformationforCaltoolsSoftware.
                                      HDUCLAS2='RSP_MATRIX', #/KeywordinformationforCaltoolsSoftware.
                                      HDUVERS='1.3.0', #/KeywordinformationforCaltoolsSoftware.
                                      HDUCLAS3='DETECTOR', #/KeywordinformationforCaltoolsSoftware.
                                  )),
            )
    
    @property
    def ebounds_hdu(self):
        return fits.BinTableHDU.from_columns([
                    fits.Column(name='CHANNEL', array=np.arange(self._e_min.shape[0]), format='1I'),
                    fits.Column(name='E_MIN', array=self._e_min, format='1E'),
                    fits.Column(name='E_MAX', array=self._e_max, format='1E'),
                ],
                header=fits.Header(cards=dict(
                                      EXTNAME='EBOUNDS', #           / name of this binary table extension
                                      TLMIN1=0, #/ First legal channel number
                                      TLMAX1=len(self._e_min), #                  511 / Highest legal channel number
                                      TELESCOP=self._telescop, #          / mission/satellite name
                                      INSTRUME=self._instrume, #           / instrument/detector
                                      #FILTER  = 'NONE    '           / filter information
                                      CHANTYPE='PI', #           / Type of channels (PHA, PI etc)
                                      DETCHANS=len(self._e_min), # / Total number of detector PHA channels
                                      #SMOOTHED=                    0 / 0 = raw, 1-12 = smooth, -1 = ep-lin, -2 = mean-
                                      HDUCLASS='OGIP', # / Keyword information for Caltools Software.
                                      HDUCLAS1='RESPONSE', # / Keyword information for Caltools Software.
                                      HDUCLAS2='EBOUNDS', # / Keyword information for Caltools Software.
                                      HDUVERS ='1.2.0', # / Keyword information for Caltools Software.
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
