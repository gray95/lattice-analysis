#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
Edited by Gaurav Ray March 2021
"""

import matplotlib
matplotlib.use('Agg')
import os
import sys
sys.path.append('./plotters')
import lsqfit
from corrfitter import Corr2,Corr3,CorrFitter
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import math as m
from math import exp as num_exp
import matplotlib.pyplot as plt
import g2tools as g2
import lsqfit as lsq
from commonweal import w0, w0overa, ZV, ZVqed
from fitting import build_models, build_prior, fitargs, make_data
from fitting import print_results

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
#-------------
a_str = 'vc' #|
#-------------
hbarc = 0.197326968
a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

ZV = ZV[a_str]
ZVqed = ZVqed[a_str]*ZV

def main(tstr):
    dfile = '../data/qqed/vcoarse/ms_vcoarse.gpl'
   
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    cdata = madedata[0]
    T = cdata.size / len(madedata[1]) 		# extent in time dir
    tag01 = madedata[1][0]
    tag02 = madedata[1][1]

    print("time extent = %d"%(T))
    print(dfile)
    print(len(cdata[tag01]))
    print("t* = %d = %sfm"%(tstr, a*hbarc*tstr))

    suggestedsvdcut = madedata[2]
    tmin = 2
    svdcut = suggestedsvdcut

    tdata = range(T)
    global tfit
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    
    pfile = None #"vector_fit.p" # last fit
    global ccdata
    ccdata = cdata[tag01]
    
    fitter = CorrFitter(models=build_models(tag01,tag02,tmin,T))

#   z0 = { 'p0' : 0.5, 
# 				 'nexp' : 5  }
#   global model
#   model_noqed, model_qed = build_models(tag01,tag02,tmin,T)
#   model = model_noqed
#   fit, z = lsq.empbayes_fit(z0, fitargs)
#    print(fit.format(True))
#   print("optimal priors WITHOUT QED: ", z)

#   mm = int( m.floor(z['nexp']) )
#   NEXP = range( mm, mm+2 ) 

#   ccdata = cdata[tag02]
#   model = model_qed
#   fit, z = lsq.empbayes_fit(z0, fitargs)
#    print(fit.format(True))
#   print("optimal priors WITH QED: ", z)

    NEXP = range(3,6)

    for nexp in NEXP:
        fit = fitter.lsqfit(data=cdata,prior=build_prior(nexp, 1.0, a.mean), p0=pfile ,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)
    tstar = tstr 

    p = fit.p
    E = np.cumsum(p['dE:vec:s'])
    amp = p['a1:vec:s']
    Eo = np.cumsum(p['dEo:vec:s'])
    ampo = p['ao:vec:s']
    Eqed =  np.cumsum(p['dE:vec:qed:s'])
    ampqed = p['a1:vec:qed:s']

    print_results(fit, ainv, ZV, ZVqed)

    tcut = 200

    newdata = {}

    tags = ['nocharge', 'charge']
    newdata[tags[0]] = cdata[tag01][:tstar+1]
    newdata[tags[1]] = cdata[tag02][:tstar+1]

    # do replacement of data with fit

    fitdatauncharged = Corr2(datatag=tags[0],tdata=range(T),a=('a1:vec:s','ao:vec:s'),b=('a1:vec:s','ao:vec:s'),dE=('dE:vec:s','dEo:vec:s'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(T):
            if index >= tcut:
                newdata[tags[0]] = np.append(newdata[tags[0]],0.)
            elif index > tstar:
                newdata[tags[0]] = np.append(newdata[tags[0]],fitdatauncharged[index])

    fitdatachargeddown = Corr2(datatag=tags[1],tdata=range(T),a=('a1:vec:qed:s','ao:vec:qed:s'),b=('a1:vec:qed:s','ao:vec:qed:s'),dE=('dE:vec:qed:s','dEo:vec:qed:s'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(T):
            if index >= tcut:
                newdata[tags[1]] = np.append(newdata[tags[1]],0.)
            elif index > tstar:
                newdata[tags[1]] = np.append(newdata[tags[1]],fitdatachargeddown[index])

#    print('THIS IS THE INFORMATION CENTRE')
#    print(len(newdata[tags[0]]), len(newdata[tags[1]]))
    print('tmin, tstar', tmin, tstar)

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamus = g2.a_mu(vpol,1/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamus = g2.a_mu(vpol,1/3.)

    amus_rt = chargedamus/unchargedamus
    amus_diff = chargedamus-unchargedamus
    amus_diff_rt = amus_diff/unchargedamus

    gv.dump(amus_rt, 'test.p')

    ## PRINT RESULTS FOR anomaly
    print('\n[s] a_mu[QCD+QED] || ratio || diff || diff_rt\n{0:20}{1:20}{2:20}{3:20}'.format(chargedamus, amus_rt, amus_diff, amus_diff_rt))    
    inputs = { 'a':a, 'ZV':ZV, 'ZVqed':ZVqed, 'a0':amp[0], 'a1':amp[1], 'm0':E[0] }
    outputs = { 'amu':chargedamus, 'amu diff':amus_diff, 'amu rt':amus_rt }
    print(gv.fmt_errorbudget(outputs=outputs, inputs=inputs))
######################################################################################
    

if __name__ == '__main__':
    for i in [17]:
        main(i)
