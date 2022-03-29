from __future__ import print_function   # makes this work for python2 and 3

import matplotlib.pyplot as plt
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import sys

from pseudoscalar_fit_params import TFIT, TDATA, NEXP, TP, savefile
from pseudoscalar_fit_params import ainv_vc, ainv_c, ainv_f
from pseudoscalar_fit_params import tags, corrpath, NMASSES

ainv = ainv_f

def fit_data(filename_in, keys, otherkeys):
    dset = cf.read_dataset(filename_in) # read data 
	
    # If we don't have many samples this next part suggests an svdcut
    s = gv.dataset.svd_diagnosis(dset, models=make_models(keys, otherkeys))
    print('svdcut =', s.svdcut) # suggested svdcut
    #s.plot_ratio(show=True)

    data = make_data(filename_in)
	
    fitter = cf.CorrFitter(models=make_models(keys, otherkeys))
    p0 = None
    for N in NEXP:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, NMASSES) 
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=s.svdcut) # ncg 
        p0 = fit.pmean
        #print(fit.format(pstyle=None if N < 10 else 'm'))

    print(fit.format(pstyle='m'))
    if fit.chi2/fit.dof > 1:
        print("chi2/dof = %f"%(fit.chi2/fit.dof))
        sys.exit(0)
    print_results(fit, NMASSES)
    # Add noise to fit to check svd cut isn't too severe
    noisyfit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=s.svdcut, add_svdnoise=True, add_priornoise=True)
    red_chi2 = noisyfit.chi2/noisyfit.dof
    print("chi2 with noise = %f"%red_chi2)
    if red_chi2 > 1.2:
        print("fit not stable enough under addition of noise")
        #sys.exit(0)

    p = fit.p
    obs = { "E:n":[] , "E:u":[], "E:d":[], "a:n":[], "a:u":[], "a:d":[] }

    for m in range(NMASSES):
        m = 'm'+str(m)
        for Q in ['n', 'u', 'd']:
            obs["E:"+Q].append(p[m+'dE:'+Q][0])
            obs["a:"+Q].append(p[m+'a:'+Q][0])

    gv.dump(obs, savefile, add_dependencies=True)
    
###### END OF MAIN() ####################
#######################################################################################

def make_data(filename):
    """ Read data, compute averages/covariance matrix for G(t). """
    return gv.dataset.avg_data(cf.read_dataset(filename))

def make_models(keys, otherkeys):
    """ Create corrfitter model for G(t). """
    models = []
    count = 0
    for k in keys['Q0']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=c+'a:n', b=c+'a:n', dE=c+'dE:n'))
        count = count + 1
    count = 0
    for k in keys['Q101']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=c+'a:d', b=c+'a:d', dE=c+'dE:d'))
        count = count + 1
    count = 0
    for k in keys['Q202']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=c+'a:u', b=c+'a:u', dE=c+'dE:u'))
        count = count + 1
    return models    

def make_prior(N, NMASSES):
    """ Create prior for N-state fit. """
    prior = collections.OrderedDict()   
    for m in range(NMASSES):
        m = 'm'+str(m)
        prior['log('+m+'a:n)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'a:u)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'a:d)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'dE:n)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dE:u)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dE:d)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
    return prior

def print_results(fit, NMASSES):
    p = fit.p
    g = collections.OrderedDict()

    for i in range(NMASSES):
        m = 'm'+str(i)
        for Q in ['n', 'u', 'd']:
            g[m+'E:'+Q] = np.cumsum(p[m+'dE:'+Q])[0]
            g[m+'a:'+Q] = p[m+'a:'+Q][0]
   
    print(gv.tabulate(g, ncol=2))

if __name__ == '__main__':
    fit_data(corrpath, tags, otherkeys=None)

