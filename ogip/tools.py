from typing import Callable, Optional, Tuple
import numpy as np

import scipy
import scipy.stats
import scipy.optimize

from ogip.spec import PHAI, RMF, ARF

def convolve(source_model: Callable, rmf: RMF, arf: Optional[ARF]=None):
    if arf is None:
        _arf = np.ones_like(rmf.c_energ)
    else:
        _arf = arf._arf    

    return (np.outer(_arf * source_model(rmf.c_energ)*rmf.d_energ,np.ones_like(rmf._e_min))*rmf._matrix).sum(0)


def fit(source_model_generator: Callable, p0: list, pha: PHAI, rmf: RMF, arf: Optional[ARF]=None, mask=None, method='COBYLA'):
    r = scipy.optimize.minimize(
        lambda p: get_mloglike(pha, source_model_generator(p), rmf=rmf, arf=arf, mask=mask), 
        p0,
        method=method)
    return r, source_model_generator(r.x)

# def get_unfolding_factor(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF]=None) -> Tuple[np.ndarray, np.ndarray]:    
def get_unfolding_factor(model: Callable, rmf: RMF, arf: Optional[ARF]=None) -> Tuple[np.ndarray, np.ndarray]:    
    model_spec = convolve(model, rmf, arf)

    e1 = rmf._e_min
    e2 = rmf._e_max
    
    return model(e1) / model_spec * rmf.d_e


def get_loglike(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF]=None, mask=None):
    convolved = convolve(model, rmf, arf)

    if mask is None:
        mask = np.ones_like(convolved, dtype=bool)
    
    return np.sum(scipy.stats.norm(convolved, pha._stat_err).logpdf(pha._rate)[mask])


def get_mloglike(*args, **kwargs):
    return -1 * get_loglike(*args, **kwargs)


def synthesise(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF]=None):
    return scipy.stats.norm(convolve(model, rmf, arf), pha._stat_err).rvs()


def plot(pha: PHAI, model: Callable, rmf: RMF, arf: Optional[ARF]=None, 
         fig=None, label_prefix=None, plot_kwargs={}, unfolded=False, e_power=0, step=False):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pylab as plt

    ie1 = rmf._energ_lo
    ie2 = rmf._energ_hi
    e1 = rmf._e_min
    e2 = rmf._e_max
    
    model_spec = convolve(model, rmf, arf)

    def prefix(l):
        if label_prefix is None:
            return l
        else:
            return f"{label_prefix}: {l}"

    if fig is None:
        plt.figure(figsize=(15,10))


    if unfolded:        
        plt.plot(rmf.c_energ, model(rmf.c_energ) * rmf.c_energ**e_power)
    else:
        plt.step(
            e1,
            model_spec/rmf.d_e,
            where='post',
            label=prefix("model"),
            **plot_kwargs
        )
    
    if pha is not None:
        if unfolded:
            unfolding_factor = get_unfolding_factor(model, rmf, arf)
        else:
            unfolding_factor = 1

        if step:
            x = plt.step(
                rmf.c_e,
                pha._rate/rmf.d_e * unfolding_factor * rmf.c_e**e_power,
                where='mid',
                label=prefix("data"),
                **plot_kwargs
            )
            color = x[0].get_color()
        else:
            color = None
                    
        plt.errorbar(
            rmf.c_e,
            pha._rate/rmf.d_e * unfolding_factor * rmf.c_e**e_power,
            pha._stat_err/rmf.d_e * unfolding_factor * rmf.c_e**e_power,
            xerr=rmf.d_e/2,
            ls="",
            c=color,
            **plot_kwargs
        )

    plt.legend()     
    plt.loglog()
            
    plt.xlim([15, 600])
    plt.grid()


def transform_rmf(rmf: RMF, arf: ARF, bias_function: Callable, preserved_projection=None, tolerance=0.05) -> RMF:
    new_rmf = RMF.from_arrays(
        energ_lo=bias_function(rmf._energ_lo.copy()),
        energ_hi=bias_function(rmf._energ_hi.copy()),
        e_min=rmf._e_min.copy(),
        e_max=rmf._e_max.copy(),
        matrix=rmf._matrix.copy()
    )

    if preserved_projection is not None:
        # projection on model should be preserved
        corr = convolve(preserved_projection, new_rmf, arf) / convolve(preserved_projection, rmf, arf)

        new_rmf._matrix /= np.outer(np.ones_like(new_rmf._energ_lo), corr)

        s = convolve(preserved_projection, rmf, arf) 
        t_s = convolve(preserved_projection, new_rmf, arf)

        assert np.all(np.abs((s - t_s)/s) < tolerance)
    

    return new_rmf

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

def le_bias_function(energy, bias_keV):
    # TODO: preserves 60 and 511, shifts LE
    return energy + bias_keV