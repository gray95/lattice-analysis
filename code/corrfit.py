from __future__ import print_function   # makes this work for python2 and 3

import matplotlib.pyplot as plt
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import sys

from parameters import TFIT, TDATA, NEXP, TP
from parameters import OSC, s_coeff

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
        print_results(fit)
        print(fit.format(pstyle=None if N < 10 else 'm'))
    

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
    prior['a'] = gv.gvar(['1.00(0.99)'] + (N-1)*['0.01(0.99)'])
    prior['log(dE)'] = gv.log(gv.gvar(['2.0(1.0)'] + (N-1)*['0.3(2)']))
    
    #----------OSC PARAMS-------------#
    if OSC:
      prior['ao'] = gv.gvar(['1.00(0.99)'] + (N-1)*['0.01(0.99)'])    
      prior['log(dEo)'] = gv.log(gv.gvar(['2.0(1.0)'] + (N-1)*['0.3(2)']))
    return prior

def print_results(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    
    if OSC:
      Eo = np.cumsum(p['dEo'])
      ao = p['ao']
	
    print('\n\nE:', end='')
    for j in range(E.shape[0]):
        print('\t', E[j], end='\t') 
        
    print('\n\na:', end='')
    for j in range(E.shape[0]):
       print('\t', a[j], end='\t')
     
    if OSC:    
      print('\n\nEo:', end='')
      for j in range(Eo.shape[0]):
          print('\t', Eo[j], end='\t')  
     
      print('\n\nao:', end='')
      for j in range(Eo.shape[0]):
          print('\t', ao[j], end='\t')
        
    print('\n=====================================================================================\n')
