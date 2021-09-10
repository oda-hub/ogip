import logging
import os
import numpy as np
import pytest

import ogip.spec

logging.basicConfig(level=logging.DEBUG)

# Roques & Jourdain 2018
def crab_ph_cm2_s_kev(en):
    K=7.417e-4
    al1=-1.98
    al2=-2.33
    Ec=500.
    f=K*(en/100)**al1*(np.exp(-en/Ec))
    m=en>Ec*(al1-al2)
    f[m]=(K*((al1-al2)*Ec/100)**(al1-al2)*(en/100)**al2*np.exp(-(al1-al2)))[m]

    return f

def test_crab(crab_spectra):
    crab_pha_fn, crab_rmf_fn, crab_arf_fn = crab_spectra
    
    crab_pha = ogip.spec.PHAI.from_file_name(crab_pha_fn)
    crab_rmf = ogip.spec.RMF.from_file_name(crab_rmf_fn)
    crab_arf = ogip.spec.ARF.from_file_name(crab_arf_fn)
    
    ie1 = crab_rmf._energ_lo
    ie2 = crab_rmf._energ_hi
    e1 = crab_rmf._e_min
    e2 = crab_rmf._e_max
    
    source=crab_ph_cm2_s_kev(ie1)
        
    model_spec = np.outer(crab_arf._arf * source*(ie2-ie1),np.ones_like(e1))*crab_rmf._matrix   
    
   
    #csource=np.outer(arf['SPECRESP']*source*(ie2-ie1),np.ones_like(rmf_eb['E_MIN']))*rmf_mt['MATRIX']
        
    
    if True:
        import matplotlib.pylab as plt
        plt.figure()
        plt.step(
            e1,
            model_spec.sum(0)/(e2 - e1),
            where='pre',
            label="model"
        )
        
        plt.step(
            e1,
            crab_pha._rate/(e2 - e1),
            where='pre',
            label="data"
        )
        
        plt.loglog()
        
        # extra_ticks(
        #     [20,25,30,35]
        # )
        
        plt.xlim([15, 600])
        plt.grid()
        plt.savefig("conv.png")

    # rate_n = np.nansum(spec.data['RATE'][(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)])

    # n = csource[:,(rmf_eb['E_MIN']>e1) & (rmf_eb['E_MAX']<e2)].sum()
    
    # erg_cm2_per_count = flux_erg_cm2_s_eband/n
    
    # print("response norm in", e1,e2,"is",n, "rate norm", rate_n, "erg/cm2 per count", erg_cm2_per_count)

    # return erg_cm2_per_count

    
