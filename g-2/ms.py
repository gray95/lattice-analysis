#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
Edited by Gaurav Ray March 2021
"""

import matplotlib
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
from commonweal import w0, w0overa, ZV, ZVqed_d, hbarc
from fitting import build_models, build_prior, fitargs, make_data
from fitting import print_results

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
#-------------
a_str = 'f' #|
#-------------
a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

ZV = ZV[a_str]
ZVqed = ZVqed_d[a_str]*ZV

base = 'fine/'
name = 'ms_rho_finep'
datafile = os.path.join('../data/qqed', base, name+'.gpl')

tmin = 2
T_STAR = [29]
bnSze = 1

store_fit = './fits/fine/ms.p'
store_res = None #'./store/fine/ms.p'

## t* reconstruction
store_corr = './fits/tstar_fine.p'

def main(tstr):
    dfile = datafile
   
    madedata = make_data(dfile,norm=3.,binsize=bnSze) # factor of 3 for colour (missed in extraction)
    data = gv.dataset.avg_data(madedata[0])
    tag01 = list(madedata[1])[0]
    tag02 = list(madedata[1])[1]
    T = data.size / len(madedata[1]) 		# extent in time dir
    T = int(T)
    print("time extent = %d"%(T))
    print("t* = %d = %sfm"%(tstr, a*hbarc*tstr))

    ### svd diagnosis ###
    s = gv.dataset.svd_diagnosis(madedata[0], models=build_models(tag01,tag02,tmin,T))
    s.plot_ratio(show=True)
    suggestedsvdcut = s.svdcut
    svdcut = suggestedsvdcut
    ####################

    print('svd cut = %f'%svdcut)

    pfile = store_fit # last fit
    
    fitter = CorrFitter(models=build_models(tag01,tag02,tmin,T))

    NEXP = range(3,6)

    for nexp in NEXP:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp, 1.7, a.mean), p0=pfile ,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit.format(pstyle='m'))
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
    newdata[tags[0]] = data[tag01][:tstar+1]
    newdata[tags[1]] = data[tag02][:tstar+1]

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
    print("amu is %s"%unchargedamus)
    print("amu with qed is %s"%chargedamus)
    print("diff is %s"%amus_diff)
    info = {"diff":amus_diff, "ratio":amus_rt, "amus":unchargedamus, "amus_qed":chargedamus}
    corr = {"corr":data, "newcorr":newdata}

    gv.dump( info, store_res, True)
    gv.dump( corr, store_corr, True)

    ## PRINT RESULTS FOR anomaly
    inputs = { 'fit':fit.p, 'a':a, 'ZV':ZV, 'ZVqed':ZVqed, 'a0':amp[0], 'a1':amp[1], 'm0':E[0] }
    outputs = { 'amu':chargedamus, 'amu diff':amus_diff, 'amu rt':amus_rt }
    #print(gv.fmt_errorbudget(outputs=outputs, inputs=inputs))
######################################################################################
    

if __name__ == '__main__':
    for i in T_STAR:
        main(i)
