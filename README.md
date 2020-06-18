# Python library to read, write, and validate "OGIP" files

* aims to support OGIP formats, especially [OGIP/92-007 (Spectral File Format)](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node5.html) ([pdf](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007.pdf))
* does not aim to implement the strictest form of the standard, instead trying to understand the content with a well-defined process
* relies only on [astropy](https://www.astropy.org/)

Existing alternatives:

* [heasoft](https://heasarc.gsfc.nasa.gov/lheasoft/) contains pyxpec and heasp modules, which read ogip files. 
* [sherpa](https://cxc.cfa.harvard.edu/sherpa/) and [3ml](https://github.com/threeML/threeML/) can read and write ogip files

The alertnatives share common disadvantage: they are complex and not easily portable packages.

## Disadvantages:

* only python, hence slow
