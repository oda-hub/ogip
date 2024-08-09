import logging
import numpy as np # type: ignore
import astropy.io.fits as fits # type: ignore

logger = logging.getLogger()


class LightCurve:

    @staticmethod
    def from_file_name(fn):
        raise NotImplementedError

    def to_long_string(self):
        return repr(self)

    ###



double_time_definition_keywords = [ ('MJDREFI', 'MJDREFF'), 
                                    ('TSTARTI', 'TSTARTF'),
                                    ('TSTOPI', 'TSTOPF'),
                                    ('TIMEZERI', 'TIMEZERF')]

single_time_definition_keywords = ['MJDREF',
                                    'TSTART',
                                    'TSTOP',
                                    'TIMEZERO',
                                    'TIMESYS', 
                                    'TIMEUNIT',
                                   'CLOCKCOR',
                                   'TIMEREF',
                                   'TIMEDEL']

optional_keywords = ['DATE-OBS' , 'DATE-END',
                                      'TELESCOP', 
                                      'INSTRUME', 
                                      'FILTER',
                                      'OBJECT',
                                      'RA',
                                      'DEC',
                                      'EQUINOX',
                                      'RADECSYS',
                                      'ORIGIN',
                                      'DATE',
                                      'AUTHOR']


class rate(LightCurve):
    _time = None
    _timedel = None
    _rate1 = None
    _err1 = None
    _fracexp = None
    _keywords = {}

    _filename = None
    
    def __init__(self, filename):

        self._keywords = {kk : None for kk in single_time_definition_keywords}

        methods_tried= {}
        for n, m in self.__dict__.items():
            if n.startswith('from_file_name_'):
                try:
                    m.__func__(filename)
                except Exception as e:
                    logger.debug("failed to read with %s.%s %s: %s", self.__name__, n, m, e)
                    methods_tried[n] = e

        raise Exception(f"failed to read with any method, tried: {methods_tried}")

    
    def read_keywords(self, hdu):
         
        for kk in single_time_definition_keywords + optional_keywords:
            if kk in hdu.header:
                self._keywords[kk] = hdu.header[kk]
                logger.debug(f'Updating {kk}')
        for cc in double_time_definition_keywords:
            if cc[0] in hdu.header and cc[1] in hdu.header:
                tmp_val = hdu.header[cc[0]] + hdu.header[cc[1]]
                for kk in single_time_definition_keywords[0:4]:
                    if cc[0].startswith(kk[0:-1]):
                        logger.debug(f'Updating {kk} from {cc[0]} and {cc[1]}')
                        self._keywords[kk] = tmp_val
        
        
    
    def from_file_name_osa(self, fn):
        f = fits.open(fn)
        issue = True
        for e in f:
            try:
                    
                    self._filename = fn
                    self._time=e.data['TIME']
                    self._timedel=e.data['TIMEDEL']
                    self._rate1=e.data['RATE']
                    self._err1=e.data['ERROR']

                    self.read_keywords(e)
                    # set if I find an extension
                    issue = False
                                    
            except Exception as ex:
                logger.debug("failed to read from %s: %s", e, ex)
        if issue:
            raise Exception("unable to read from any extension")

    def from_file_name_normal(self, fn):
        f = fits.open(fn)
        try:
            hdu = f['RATE']
        except:
            logger.error('Unable to find the extenision RATE in ' + fn)
            raise Exception('Unable to find the extenision RATE in ' + fn)
        
        self.read_keywords(hdu)
        
        if 'TIME' in hdu.data.names:
            self._time = hdu.data['TIME']
        else:
            logger.warning("No time column, use keywords")
            self._time = np.linspace(self._keywords['TSTART'],self._keywords['TSTOP'], 
                                      self._keywords['TIMEDEL']  )
        
        #char read_dx = 0; //check if need to read TIMEDEL, 0 no read dx, 1 read TIMEDEL keyword, 2 read XAX_E (ignoring TIMEDEL keyword
        if 'XAX_E' in hdu.data.names:
            self._timedel = hdu.data['XAX_E'] * 2
        elif 'TIMEDEL' in hdu.data.names:
            self._timedel = hdu.data['TIMEDEL']
        else:
            self._timedel = self._time[1:] - self._time[0:-1]
            self._timedel = np.append(self._timedel, self._timedel[-1])
        
        if 'RATE' in hdu.data.names:
            self._rate1 = hdu.data['RATE']
        elif 'COUNTS' in hdu.data.names:
            self._rate1 = hdu.data['COUNTS']
        else:
            rate_columns = []
            for kk in hdu.data.names:
                if kk.startswith('RATE'):
                    rate_columns.append(kk)

            if len(rate_columns) == 0 :
                logger.error(f'No RATE or COUNTS in {fn}')
                raise Exception(f'No RATE or COUNTS in {fn}')
            else:
                for kk in rate_columns:
                    self.__setattr__(f'_{kk.lower()}', hdu.data[kk])
        
        if 'ERROR' in hdu.data.names:
            self._err1 = hdu.data['ERROR']
        else:
            err_columns = []
            for kk in hdu.data.names:
                if kk.startswith('ERROR'):
                    err_columns.append(kk)

            if len(err_columns) == 0 :
                logger.warning(f'No Errors in {fn} assuming Poissonian')
                for kk_r in rate_columns:
                    err_kk = kk_r.replace('RATE', 'ERROR')
                    self.__setattr__(f'_{err_kk.lower()}', np.sqrt(hdu.data[kk_r]))
                
            else:
                for kk in err_columns:
                    self.__setattr__(f'_{kk.lower()}', hdu.data[kk])
            
        
        if 'FRACEXP' in hdu.data.names:
            self._fracexp = hdu.data['FRACEXP']
        else:
            self._fracexp = np.ones(len(self._rate))


    @property
    def lc_hdu(self):
        # https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/rates/ogip_93_003/ogip_93_003.html

        columns = [
                    fits.Column(name='RATE', array=self._rate, format='1E'),
                    fits.Column(name='ERROR', array=self._stat_err, format='1E'),
                    fits.Column(name='TIME', array=self._sys_err, format='D'),
                    fits.Column(name='TIMEDEL', array=self._timedel, format='D'),
                    fits.Column(name='FRACEXP', array=self._fracexp, format='1E'),
                ]
        header = fits.Header(cards=dict(
                                      EXTNAME="RATE", # - the name (i.e. type) of the extension
                                      TELESCOP="", # - the "telescope" (i.e. mission/satellite name).
                                      INSTRUME="", #- the instrument/detector.
                                      FILTER="", #- the instrument filter in use (if any)
                                      HDUCLASS="OGIP", # - should contain the string "OGIP" to indicate that this is an OGIP style file.
                                      HDUCLAS1="LIGHT CURVE", # -
                                      HDUCLAS2="RATE", # -
                                      TIMVERSN= 'OGIP/93-003', #OGIP memo number where the convention used
                                      OBJECT  = '', #   / name of the observed object
                                      RA      = '', #                  / source right ascension in degree 
                                      DEC     = '', #                   / source declination in degree 
            # RA--NOM =                    / r.a. nominal pointing in degree
            # DEC--NOM=                    / dec. nominal pointing in degree
                                      EQUINOX = 2000.0, #             / equinox for ra and dec
                                      RADECSYS= 'FK5',#              / world coord. system for this file (FK5 or FK4)
                                      ORIGIN  = '',#         / who produced the fits file
                                      DATE    = '', #         / when the fits file was created
                                      AUTHOR  = '',   #  /name of program that produced this file
                                      ))
        for kk in ['DATE-OBS' , 'DATE-END']:
            header[kk] = ''
        for kk, vv in self._keywords:
            header[kk] = vv

        return fits.BinTableHDU.from_columns(columns, header=header)

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


