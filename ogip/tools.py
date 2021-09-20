import numpy as np
from ogip.spec import PHAI, RMF, ARF

def convolve(source_model: callable, rmf: RMF, arf: ARF=None):
    ie1 = rmf._energ_lo
    ie2 = rmf._energ_hi
    e1 = rmf._e_min
    e2 = rmf._e_max
    
    if arf is None:
        _arf = np.ones_like(ie1)
    else:
        _arf = arf._arf    

    return (np.outer(_arf * source_model(ie1)*(ie2-ie1),np.ones_like(e1))*rmf._matrix).sum(0)

def plot(pha: PHAI, model: callable, rmf: RMF, arf: ARF=None, fig=None, label_prefix=None, plot_kwargs={}):
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

    plt.step(
        e1,
        model_spec/(e2 - e1),
        where='post',
        label=prefix("model"),
        **plot_kwargs
    )
    
    if pha is not None:
        x = plt.step(
            e1,
            pha._rate/(e2 - e1),
            where='post',
            label=prefix("data"),
            **plot_kwargs
        )
                
        plt.errorbar(
            (e1 * e2)**0.5,
            pha._rate/(e2 - e1),
            pha._stat_err/(e2 - e1),
            ls="",
            c=x[0].get_color(),
            **plot_kwargs
        )


    plt.legend()     
    plt.loglog()
    
        
    plt.xlim([15, 600])
    plt.grid()

def transform_rmf(rmf: RMF, arf: ARF, bias_function: callable, preserved_projection=None, tolerance=0.05) -> RMF:
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