from __future__ import print_function   # makes this work for python2 and 3

import matplotlib.pyplot as plt
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import sys

from parameters import TFIT, TDATA, NEXP, TP
from parameters import OSC, s_coeff, ainv

def fit_data(filename_in, key, otherkey):
    dset = cf.read_dataset(filename_in) # read data 
	
	# If we don't have many samples this next part suggests an svdcut
    #s = gv.dataset.svd_diagnosis(dset, models=make_models(key, otherkey, n))
    #print('svdcut =', s.svdcut) # suggested svdcut
    #s.plot_ratio(show=True)
	############################################################################################
	
    data = make_data(filename_in)
	
    fitter = cf.CorrFitter(models=make_models(key, otherkey))
    p0 = None
    for N in NEXP:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N) 
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0) # add_svdnoise=True, add_priornoise=True, svdcut=s.svdcut
        p0 = fit.pmean
        print(fit.format(pstyle=None if N < 10 else 'm'))
    print_results(fit)

def make_data(filename):
    """ Read data, compute averages/covariance matrix for G(t). """
    return gv.dataset.avg_data(cf.read_dataset(filename))

def make_models(key, otherkey):
    """ Create corrfitter model for G(t). """
    if OSC:
      return [cf.Corr2(datatag=key, otherdata=otherkey, tdata=TDATA,tfit=TFIT, tp=TP, a=('a','ao'),b=('a','ao'),dE=('dE','dEo'), s=s_coeff)]
    else:
      return [cf.Corr2(datatag=key, otherdata=otherkey, tdata=TDATA, tfit=TFIT, tp=TP, a='a', b='a', dE='dE')]


def make_prior(N):
    """ Create prior for N-state fit. """
    prior = collections.OrderedDict()    
    prior['log(a)'] = gv.log(gv.gvar(['1.00(0.99)'] + (N-1)*['0.01(0.99)']))
    prior['log(dE)'] = gv.log(gv.gvar(['1.4(1), 0.1(1)'] + (N-2)*['0.5(5)']))
    
    #----------OSC PARAMS-------------#
    if OSC:
      prior['log(ao)'] = gv.log(gv.gvar(['1.00(0.99)'] + (N-1)*['0.01(0.99)']))    
      prior['log(dEo)'] = gv.log(gv.gvar(['2.0(5)'] + (N-1)*['0.5(5)']))
    return prior

def print_results(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    
    if OSC:
      Eo = np.cumsum(p['dEo'])
      ao = p['ao']
	
    print('\n\n\tE (GeV)\t\t\ta')
    print('---------------------------------------------------------')
    for j in range(E.shape[0]):
        print(" %d:    %s \t\t%s" % (j, ainv*E[j], a[j])) 
        
    if OSC:    
        print('\n\n\tEo (GeV)\t\t\tao')
        print('---------------------------------------------------------')
        for j in range(Eo.shape[0]):
            print(" %d:    %s \t\t%s" % (j, ainv*Eo[j], ao[j])) 
       
    print('\n=====================================================================================\n')
