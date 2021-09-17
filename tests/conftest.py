import os
import pytest

@pytest.fixture
def crab_spectra():
    crab_pha_fn = 'crab_pha.fits'
    crab_rmf_fn = 'crab_rmf.fits'
    crab_arf_fn = 'crab_arf.fits'
    
    if os.path.exists(crab_pha_fn) and os.path.exists(crab_rmf_fn) and os.path.exists(crab_arf_fn):
        pass
    else:
        import oda_api.api
        disp = oda_api.api.DispatcherAPI()
        prod = disp.get_product(
            instrument='isgri',
            product='isgri_spectrum',
            scw_list="066500220010.001",
            osa_version="OSA10.2"
            )

        
        prod.isgri_spectrum_0_Crab_isgri_spectrum.write_fits_file(crab_pha_fn)
        prod.isgri_spectrum_1_Crab_isgri_arf.write_fits_file(crab_arf_fn)
        prod.isgri_spectrum_2_Crab_isgri_rmf.write_fits_file(crab_rmf_fn)

    return crab_pha_fn, crab_rmf_fn, crab_arf_fn