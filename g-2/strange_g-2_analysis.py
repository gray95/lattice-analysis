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

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
#-------------
a_str = 'c' #|
#-------------
hbarc = 0.197326968
a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing (fm): ", a*hbarc)

ZV = ZV[a_str]
ZVqed = ZVqed[a_str]*ZV

def main(tstr):
    dfile = '../data/qqed/coarse/ms_rho_coarse.gpl'
   
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    cdata = madedata[0]
    T = cdata.size / len(madedata[1]) 		# extent in time dir
    tag01 = madedata[1][0]
    tag02 = madedata[1][1]
    print("time extent is: ", T)
    print(dfile)
    print(len(cdata[tag01]))
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
    amu_diff = chargedamus-unchargedamus
    amu_diff_rt = amu_diff/unchargedamus

    ## PRINT RESULTS FOR anomaly
    print('\n[s] a_mu[QCD+QED] || ratio || diff || diff_rt\n{0:20}{1:20}{2:20}{3:20}'.format(chargedamus, amus_rt, amu_diff, amu_diff_rt))    
    #print('[s] a_mu[QCD+QED]{0:>20}\n|a_mu[QCD+QED]-a_mu[QCD]{1:>20}\n|a_mu/a_mu {2:>20}'.format(chargedamus, amu_rt, amu_diff, amu_diff_rt))
######################################################################################
    
def build_prior(nexp, dp, a):
    prior = gv.BufferDict()
    nexp = int(nexp)
    prior.add('log(a1:vec:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(a1:vec:s)'][0] = log(gv.gvar(0.5, 0.5*dp))
    prior.add('log(ao:vec:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(ao:vec:s)'][0] = log(gv.gvar(0.5, 0.5*dp))
    prior.add('log(dE:vec:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:s)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:s)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:qed:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(a1:vec:qed:s)'][0] = log(gv.gvar(0.5,0.5*dp))
    prior.add('log(ao:vec:qed:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(ao:vec:qed:s)'][0] = log(gv.gvar(0.5,0.5*dp))
    prior.add('log(dE:vec:qed:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:qed:s)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:qed:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:qed:s)'][0] = log(gv.gvar(a,a*dp))
   #print(prior)
    return prior
##

def fitargs(z):
    dp = z['p0']
    nexp = z['nexp']		
    prior = build_prior(nexp, dp)    
    return dict(prior=prior, fcn=model.fitfcn, data=ccdata[2:49])
##

def build_models(tag01,tag02,tmin,T):
 
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T # periodic
    models = [
        Corr2(datatag=tag01,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:s','ao:vec:s'),b=('a1:vec:s','ao:vec:s'),
            dE=('dE:vec:s','dEo:vec:s'),s=(1.,-1.)), 
        Corr2(datatag=tag02,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:s','ao:vec:qed:s'),b=('a1:vec:qed:s','ao:vec:qed:s'),
            dE=('dE:vec:qed:s','dEo:vec:qed:s'),s=(1.,-1.)),
    ]
    return models[0], models[1]
##    

def make_data(filename,norm=1):
    dset = gv.dataset.Dataset(filename)
    dset = gv.dataset.bin_data(dset, binsize=1)
    s = gv.dataset.svd_diagnosis(dset)
    print(dset.keys())
    for tag in dset.keys():
        dset[tag] = norm*np.array(dset[tag])
        #dset[tag] = dset[tag][:3]
        print("number of cfgs| ", tag, " : ", dset[tag].shape[0]) 
    cdata = gv.dataset.avg_data(dset)
    return (cdata, dset.keys(), s.svdcut)

def print_results(fit, ainv, Zv,Zvqed):
    p = fit.p
    E = np.cumsum(p['dE:vec:s'])
    a = p['a1:vec:s']
    Eo = np.cumsum(p['dEo:vec:s'])
    ao = p['ao:vec:s']


    Eqed =  np.cumsum(p['dE:vec:qed:s'])
    aqed = p['a1:vec:qed:s']

    print("Decay constant calculation")
    print("Ground Mass [QCD] = " , E[0], "lattice units")
    print("Ground Mass [QCD+QED] = " , Eqed[0], "lattice units")
    ddd = Eqed[0] - E[0]
    print("Mass diff = " , ddd , "lattice units")
    ddd *= 1000*ainv
    print("Mass diff = " , ddd , "MeV")

    print("Ground Mass [QCD] = " , E[0] * ainv, " GeV")
    print("Ground Mass [QCD+QED] = " , Eqed[0] * ainv, " GeV")
    print("Amplitide = " , a[0])

    f_phi_latt = a[0] * (2.0 / E[0])**(1/2)
    f_phi_gev = f_phi_latt * ainv * Zv 

    print("Decay constant [QCD] = " , f_phi_latt  , "lattice")
    print("Decay constant [QCD] = " , f_phi_gev  , "GeV")

    # ZVqed
    f_phi_qed_latt = aqed[0] * (2.0 / Eqed[0])**(1/2)
    f_phi_qed_gev = f_phi_qed_latt * ainv * Zvqed

    print("Decay constant [QCD+QED] = " , f_phi_qed_gev  , "GeV")

    f_rat = f_phi_qed_gev/f_phi_gev 
    print("Ratio Decay constant [QCD+QED]/[QCD] = " , f_rat )
    
    f_diff = f_phi_qed_gev - f_phi_gev 
    print("Decay constant [diff] = " , f_diff, "GeV" )

if __name__ == '__main__':
    for i in range(18,19):
        main(i)
