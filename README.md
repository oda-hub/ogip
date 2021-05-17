# Python library to read, write, and validate "OGIP" files

[![Build Status](https://travis-ci.com/volodymyrss/ogip.svg?branch=master)](https://travis-ci.com/volodymyrss/ogip)[![codecov](https://codecov.io/gh/volodymyrss/ogip/branch/master/graph/badge.svg)](https://codecov.io/gh/volodymyrss/ogip)



* aims to support OGIP formats:
  * [OGIP/92-007 (The OGIP Spectral File Format)](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html)
  * [OGIP Calibration Memo CAL/GEN/92-002 (The Calibration Requirements for Spectral Analysis, Definition of RMF and ARF file formats)](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html)
  * [OGIP Memo OGIP/93-003 (The Proposed Timing FITS File Format for High Energy Astrophysics Data)](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/rates/ogip_93_003/ogip_93_003.html).
* does not aim to implement the strictest form of the standard, instead trying to understand the content from the least information provided, following [Postel's](https://en.wikipedia.org/wiki/Robustness_principle) [Law](https://tools.ietf.org/html/rfc1122) (which has proven to be successful in TCP). 
* well-defined, versioned, and provenance-committed process
* **it just works, always** while it is also **always possible to recover what exactly was done and why**.
* relies only on [astropy](https://www.astropy.org/)

Existing alternatives:

* [heasoft](https://heasarc.gsfc.nasa.gov/lheasoft/) contains pyxpec and heasp modules, which read ogip files. 
* [sherpa](https://cxc.cfa.harvard.edu/sherpa/) and [3ml](https://github.com/threeML/threeML/) can read and write ogip files

The alertnatives share common disadvantage: they are complex and not easily portable packages. 
Also, they implement different variations of the standards, and we have a different aim here (see above).

## Disadvantages:

* only python, hence slow
