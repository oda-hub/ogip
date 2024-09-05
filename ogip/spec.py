import logging
import numpy as np  # type: ignore
import astropy.io.fits as fits  # type: ignore

logger = logging.getLogger()


def log_bins(n_bins_per_decade, global_e_min, global_e_max):
    n_bin_edges_per_decade = n_bins_per_decade + 1

    new_e_bins = np.logspace(
        np.log10(global_e_min),
        np.log10(global_e_max),
        int(np.log10(global_e_max / global_e_min) * n_bin_edges_per_decade),
    )

    logger.debug("new e bins: %s", new_e_bins)

    return new_e_bins


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
        methods_tried = {}
        for n, m in cls.__dict__.items():
            if n.startswith("from_file_name_"):
                try:
                    return m.__func__(fn)
                except Exception as e:
                    logger.debug(
                        "failed to read with %s.%s %s: %s", cls.__name__, n, m, e
                    )
                    methods_tried[n] = e

        raise Exception(
            f"failed to read with any method, tried: {methods_tried}")


class PHAI(Spectrum):
    _rate = None
    _sys_err = None
    _stat_err = None
    _quality = None

    filename = None
    _exposure = None

    @staticmethod
    def from_file_name(fn):
        return Readable.from_file_name(PHAI, fn)

    @staticmethod
    def from_file_name_osa(fn):
        f = fits.open(fn)
        for e in f:
            try:
                return PHAI.from_arrays(
                    exposure=e.header["EXPOSURE"],
                    rate=e.data["RATE"],
                    stat_err=e.data["STAT_ERR"],
                    sys_err=e.data["SYS_ERR"],
                    filename=fn,
                )
            except Exception as ex:
                logger.debug("failed to read from %s: %s", e, ex)

        raise Exception("unable to read from any extension in OSA")

    @staticmethod
    def from_file_name_normal(fn):
        logger.debug('Reading OGIP format')

        f = fits.open(fn)

        exposure = f['SPECTRUM'].header['EXPOSURE']
        if 'RATE' in f['SPECTRUM'].data.names:
            rate = f['SPECTRUM'].data['RATE']
        elif 'COUNTS' in f['SPECTRUM'].data.names:
            rate = f['SPECTRUM'].data['COUNTS'] / exposure
        else:
            rate = None

        if 'STAT_ERR' in f['SPECTRUM'].data.names:
            stat_err = f['SPECTRUM'].data['STAT_ERR']
        elif f['SPECTRUM'].header.get('POISSERR', False):
            logger.debug('Reading Poissonian uncertainties')
            if 'COUNTS' in f['SPECTRUM'].data.names:
                stat_err = np.sqrt(f['SPECTRUM'].data['COUNTS']) / exposure
            else:
                stat_err = None
        else:
            stat_err = None

        if 'SYS_ERR' in f['SPECTRUM'].data.names:
            sys_err = f['SPECTRUM'].data['SYS_ERR']
        else:
            sys_err = None

        if 'GROUPING' in f['SPECTRUM'].data.names:
            grouping = f['SPECTRUM'].data['GROUPING']
        else:
            grouping = None

        return PHAI.from_arrays(
            exposure=exposure,
            rate=rate,
            stat_err=stat_err,
            sys_err=sys_err,
            filename=fn,
            grouping=grouping
        )

    @staticmethod
    def from_arrays(
        exposure,
        rate=None,
        stat_err=None,
        sys_err=None,
        quality=None,
        counts=None,
        filename=None,
        grouping=None
    ):
        self = PHAI()

        self.filename = filename
        self._exposure = exposure

        if counts is not None:
            logger.warning("counts found: converting to rate")
            rate = counts / exposure

            if stat_err is None:
                stat_err = counts**0.5
                logger.warning(
                    "assuming poisson stat errors. Is it really what you want?"
                )

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

        if grouping is not None:
            self._grouping = grouping
        else:
            self._grouping = np.zeros_like(rate)

        return self

    @property
    def spectrum_hdu(self):
        return fits.BinTableHDU.from_columns(
            [
                fits.Column(name="RATE", array=self._rate, format="1E"),
                fits.Column(name="STAT_ERR",
                            array=self._stat_err, format="1E"),
                fits.Column(name="SYS_ERR", array=self._sys_err, format="1E"),
                fits.Column(name="QUALITY", array=self._quality, format="1I"),
            ],
            # https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node6.html
            header=fits.Header(
                cards=dict(
                    # - the name (i.e. type) of the extension
                    EXTNAME="SPECTRUM",
                    # - the "telescope" (i.e. mission/satellite name).
                    TELESCOP="",
                    INSTRUME="",  # - the instrument/detector.
                    # FILTER - the instrument filter in use (if any)
                    AREASCAL=1.0,
                    BACKSCAL=1.0,
                    # the integration time (in seconds) for the PHA data (assumed to be corrected for deadtime, data drop-outs etc. )
                    EXPOSURE=self._exposure,
                    # - the name of the corresponding background file (if any)
                    BACKFILE="",
                    # CORRFILE - the name of the corresponding correction file (if any)
                    # CORRSCAL - the correction scaling factor.
                    # - the name of the corresponding (default) redistribution matrix file (RMF; see George et al. 1992a).
                    RESPFILE="",
                    # - the name of the corresponding (default) ancillary response file (ARF; see George et al. 1992a).
                    ANCRFILE="",
                    # - should contain the string "OGIP" to indicate that this is an OGIP style file.
                    HDUCLASS="OGIP",
                    # - should contain the string "SPECTRUM" to indicate this is a spectrum.
                    HDUCLAS1="SPECTRUM",
                    # - the version number of the format (this document describes version 1.2.1)
                    HDUVERS="1.2.1",
                    # POISSERR #- whether Poissonian errors are appropriate to the data (see below). # noqa: E800
                    # CHANTYPE #- whether the channels used in the file have been corrected in anyway (see below). # noqa: E800
                    DETCHANS=len(
                        self._rate
                    ),  # - the total number of detector channels available.
                )
            ),
        )

    def to_fits(self, fn: str):
        logger.info("store to fits %s %s", self, fn)

        fits.HDUList(
            [
                fits.PrimaryHDU(),
                self.spectrum_hdu,
            ]
        ).writeto(fn, overwrite=True)

    @property
    def total_counts(self):
        return np.nansum(self._rate)

    def to_long_string(self):
        return f"{self.__class__.__name__}: exposure {self._exposure} s, total count/s {self.total_counts}, {len(self._rate)} channels "


class PHAII(Spectrum):
    pass


class RMF:
    _telescop = "not-a-telescope"
    _instrume = "not-an-instrument"

    _energ_lo = None  # type: ignore
    _energ_hi = None  # type: ignore
    _matrix = None  # type: ignore
    _e_min = None  # type: ignore
    _e_max = None  # type: ignore

    def to_long_string(self):
        return f"{self.__class__.__name__}: {len(self._e_min)} x {len(self._energ_lo)}  channels "

    @staticmethod
    def from_file_name(fn):
        return Readable.from_file_name(RMF, fn)

    @staticmethod
    def from_file_name_osaisgri(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
            energ_lo=f["ISGR-RMF.-RSP"].data["ENERG_LO"],
            energ_hi=f["ISGR-RMF.-RSP"].data["ENERG_HI"],
            matrix=np.stack(f["ISGR-RMF.-RSP"].data["MATRIX"]),
            e_min=f["ISGR-EBDS-MOD"].data["E_MIN"],
            e_max=f["ISGR-EBDS-MOD"].data["E_MAX"],
        )

    @staticmethod
    def from_file_name_osaspi(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
            energ_lo=f["SPI.-RMF.-RSP"].data["ENERG_LO"],
            energ_hi=f["SPI.-RMF.-RSP"].data["ENERG_HI"],
            matrix=np.stack(f["SPI.-RMF.-RSP"].data["MATRIX"]),
            e_min=f["SPI.-EBDS-SET"].data["E_MIN"],
            e_max=f["SPI.-EBDS-SET"].data["E_MAX"],
        )

    @staticmethod
    def from_file_name_normal(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
            energ_lo=f["MATRIX"].data["ENERG_LO"],
            energ_hi=f["MATRIX"].data["ENERG_HI"],
            matrix=np.vstack(f["MATRIX"].data["MATRIX"]),
            e_min=f["EBOUNDS"].data["E_MIN"],
            e_max=f["EBOUNDS"].data["E_MAX"],
        )

    @staticmethod
    def from_file_name_alt_normal(fn):
        f = fits.open(fn)

        return RMF.from_arrays(
            energ_lo=f["SPECRESP MATRIX"].data["ENERG_LO"],
            energ_hi=f["SPECRESP MATRIX"].data["ENERG_HI"],
            matrix=np.vstack(f["SPECRESP MATRIX"].data["MATRIX"]),
            e_min=f["EBOUNDS"].data["E_MIN"],
            e_max=f["EBOUNDS"].data["E_MAX"],
        )

    @staticmethod
    def from_arrays(energ_lo, energ_hi, matrix, e_min, e_max):
        self = RMF()

        if not (len(energ_lo) == len(energ_hi) == matrix.shape[0]):
            raise Exception(
                f"incompatible dimensions of mc energy, bounds {len(energ_lo)} {len(energ_hi)} but matrix {matrix.shape[0]}!"
            )

        if not (len(e_min) == len(e_max) == matrix.shape[1]):
            raise Exception(
                f"incompatible dimensions of channel energy, bounds {len(e_min)} {len(e_max)} but matrix {matrix.shape[1]}!"
            )

        self._energ_lo = energ_lo  # type: ignore
        self._energ_hi = energ_hi  # type: ignore
        self._matrix = matrix  # type: ignore
        self._e_min = e_min  # type: ignore
        self._e_max = e_max  # type: ignore

        return self

    @property
    def matrix_hdu(self):
        return fits.BinTableHDU.from_columns(
            [
                fits.Column(name="ENERG_LO",
                            array=self._energ_lo, format="1E"),
                fits.Column(name="ENERG_HI",
                            array=self._energ_hi, format="1E"),
                fits.Column(
                    name="N_GRP", array=np.ones_like(self._energ_lo), format="1I"
                ),
                fits.Column(
                    name="F_CHAN", array=0 * np.ones_like(self._energ_lo), format="1I"
                ),
                fits.Column(
                    name="N_CHAN",
                    array=len(self._e_min) * np.ones_like(self._energ_lo),
                    format="1I",
                ),
                fits.Column(name="MATRIX", array=self._matrix, format="PE"),
                # fits.Column(name='MATRIX', array=self._matrix, format=f'{len(self._e_min)}E'), # noqa: E800
            ],
            header=fits.Header(
                cards=dict(
                    # https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc7.1.1
                    EXTNAME="MATRIX",  # / name of this binary table extension
                    # BITPIX  =                    8 / 8-bit bytes
                    # NAXIS   =                    2 / 2-dimensional binary table
                    # NAXIS1  =                   34 / width of table in bytes
                    # NAXIS2  =                 1180 / number of rows in table
                    # PCOUNT  =              1031160 / Number of bytes acumulated in heap
                    # GCOUNT  =                    1 / one data group (required keyword)
                    # TFIELDS =                    6 / number of fields in each row
                    # TTYPE1  = 'ENERG_LO'           / label for field   1
                    # TFORM1  = 'E       '           / data format of field: 4-byte REAL
                    # TUNIT1  = 'keV     '           / physical unit of field
                    # TTYPE2  = 'ENERG_HI'           / label for field   2
                    # TFORM2  = 'E       '           / data format of field: 4-byte REAL
                    # TUNIT2  = 'keV     '           / physical unit of field
                    # TTYPE3  = 'N_GRP   '           / label for field   3
                    # TFORM3  = 'I       '           / data format of field: 2-byte INTEGER
                    # TTYPE4  = 'F_CHAN  '           / label for field   4
                    # TFORM4  = 'PI(2)   '           / data format of field: variable length array
                    # TTYPE5  = 'N_CHAN  '           / label for field   5
                    # TFORM5  = 'PI(2)   '           / data format of field: variable length array
                    # TTYPE6  = 'MATRIX  '           / label for field   6
                    # TFORM6  = 'PE(418) '           / data format of field: variable length array
                    TLMIN4=0,  # / First legal channel number
                    TLMAX4=len(self._e_min),  # / Highest legal channel number
                    TELESCOP=self._telescop,  # / mission/satellite name
                    INSTRUME=self._instrume,  # / instrument/detector
                    # FILTER  = 'NONE    '           / filter information
                    CHANTYPE="PI",  # / Type of channels (PHA, PI etc)
                    DETCHANS=len(
                        self._e_min
                    ),  # / Total number of detector PHA channels
                    LO_THRES=1.00e-07,  # / Lower probability density threshold for matrix
                    HDUCLASS="OGIP",  # /KeywordinformationforCaltoolsSoftware.
                    # /KeywordinformationforCaltoolsSoftware.
                    HDUCLAS1="RESPONSE",
                    # /KeywordinformationforCaltoolsSoftware.
                    HDUCLAS2="RSP_MATRIX",
                    HDUVERS="1.3.0",  # /KeywordinformationforCaltoolsSoftware.
                    # /KeywordinformationforCaltoolsSoftware.
                    HDUCLAS3="DETECTOR",
                )
            ),
        )

    @property
    def ebounds_hdu(self):
        return fits.BinTableHDU.from_columns(
            [
                fits.Column(
                    name="CHANNEL", array=np.arange(self._e_min.shape[0]), format="1I"
                ),
                fits.Column(name="E_MIN", array=self._e_min, format="1E"),
                fits.Column(name="E_MAX", array=self._e_max, format="1E"),
            ],
            header=fits.Header(
                cards=dict(
                    EXTNAME="EBOUNDS",  # / name of this binary table extension
                    TLMIN1=0,  # / First legal channel number
                    TLMAX1=len(
                        self._e_min
                    ),  # 511 / Highest legal channel number
                    TELESCOP=self._telescop,  # / mission/satellite name
                    INSTRUME=self._instrume,  # / instrument/detector
                    # FILTER  = 'NONE    '           / filter information
                    CHANTYPE="PI",  # / Type of channels (PHA, PI etc)
                    DETCHANS=len(
                        self._e_min
                    ),  # / Total number of detector PHA channels
                    # SMOOTHED=                    0 / 0 = raw, 1-12 = smooth, -1 = ep-lin, -2 = mean-
                    # / Keyword information for Caltools Software.
                    HDUCLASS="OGIP",
                    # / Keyword information for Caltools Software.
                    HDUCLAS1="RESPONSE",
                    # / Keyword information for Caltools Software.
                    HDUCLAS2="EBOUNDS",
                    # / Keyword information for Caltools Software.
                    HDUVERS="1.2.0",
                )
            ),
        )

    def to_fits(self, fn):
        logger.info("store to fits %s %s", self, fn)

        fits.HDUList(
            [
                fits.PrimaryHDU(),
                self.ebounds_hdu,
                self.matrix_hdu,
            ]
        ).writeto(fn, overwrite=True)

    @property
    def d_e(self):
        return self._e_max - self._e_min

    @property
    def d_energ(self):
        return self._energ_hi - self._energ_lo

    @property
    def c_e(self):
        return (self._e_max + self._e_min) / 2

    @property
    def c_energ(self):
        return (self._energ_hi + self._energ_lo) / 2


def rebin(pha: PHAI, rmf: RMF, new_e_bins) -> (PHAI, RMF):
    new_e_bins_assigned = []

    if pha._rate.shape[0] != rmf._matrix[0].shape[0]:
        raise RuntimeError()

    for new_bin_req_e1, new_bin_req_e2 in zip(new_e_bins[:-1], new_e_bins[1:]):
        m = rmf._e_min >= new_bin_req_e1
        m &= rmf._e_min < new_bin_req_e2  # really e_min here

        new_e_bins_assigned.append(
            dict(
                mask=m,
                matrix_row=rmf._matrix[:, m].sum(1),
                e_min=np.array(rmf._e_min)[m].min(),
                e_max=np.array(rmf._e_max)[m].max(),
                rate=pha._rate[m].sum(),
                stat_err=np.sum(pha._stat_err[m] ** 2) ** 0.5,
                sys_err=np.sum(pha._sys_err[m] ** 2) ** 0.5,
            )
        )

        logger.info(
            "for %s - %s new bin assigned: %s",
            new_bin_req_e1,
            new_bin_req_e2,
            {k: v for k,
                v in new_e_bins_assigned[-1].items() if len(str(v)) < 50},
        )

    return (
        PHAI.from_arrays(
            exposure=pha._exposure,
            rate=np.array([r["rate"] for r in new_e_bins_assigned]),
            stat_err=np.array([r["stat_err"] for r in new_e_bins_assigned]),
            sys_err=np.array([r["stat_err"] for r in new_e_bins_assigned]),
        ),
        RMF.from_arrays(
            energ_lo=rmf._energ_lo.copy(),
            energ_hi=rmf._energ_hi.copy(),
            matrix=np.vstack(
                [r["matrix_row"] for r in new_e_bins_assigned]
            ).transpose(),
            e_min=np.array([r["e_min"] for r in new_e_bins_assigned]),
            e_max=np.array([r["e_max"] for r in new_e_bins_assigned]),
        ),
    )


class ARF:
    _arf = None

    @staticmethod
    def from_file_name(fn):
        return Readable.from_file_name(ARF, fn)

    @staticmethod
    def from_file_name_normal(fn):
        f = fits.open(fn)

        return ARF.from_arrays(
            energ_lo=f["SPECRESP"].data["ENERG_LO"],
            energ_hi=f["SPECRESP"].data["ENERG_HI"],
            arf=f["SPECRESP"].data["SPECRESP"],
        )

    @staticmethod
    def from_file_name_osa(fn):
        f = fits.open(fn)

        return ARF.from_arrays(
            energ_lo=f["ISGR-ARF.-RSP"].data["ENERG_LO"],
            energ_hi=f["ISGR-ARF.-RSP"].data["ENERG_HI"],
            arf=f["ISGR-ARF.-RSP"].data["SPECRESP"],
        )

    @staticmethod
    def from_arrays(energ_lo, energ_hi, arf):
        self = ARF()

        if not (len(energ_lo) == len(energ_hi) == len(arf)):
            raise Exception(
                f"incompatible dimensions of mc energy, bounds {len(energ_lo)} {len(energ_hi)} but matrix {len(arf)}!"
            )

        self._energ_lo = energ_lo  # type: ignore
        self._energ_hi = energ_hi  # type: ignore
        self._arf = arf  # type: ignore

        return self

    def to_long_string(self):
        return (
            f"{self.__class__.__name__}: {len(self._arf)} energies {np.max(self._arf)}"
        )
