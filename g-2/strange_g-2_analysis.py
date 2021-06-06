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
import lsqfit
from corrfitter import Corr2,Corr3,CorrFitter
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import math
from math import exp as num_exp
import matplotlib.pyplot as plt
import g2tools as g2

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'

w0 = gv.gvar('0.1715(9)')

#w0overa = gv.gvar('1.1367(5)')	# very coarse ensemble
#ZV = gv.gvar('0.9837(20)') # vc

w0overa = gv.gvar('1.4149(6)')	# coarse ensemble
ZV = gv.gvar('0.99220(40)') # coarse - at strange mass

ZVqed = gv.gvar('0.999544(14)')*ZV # vc/c? - where are these listed

hbarc = 0.197326968
a = (w0/w0overa)/hbarc		# in units of (GeV)^-1

print("lattice spacing: ", (w0/w0overa))

def main(tstr):
    dfile = '/home/gray/Desktop/lattice-analysis/data/qed/coarse/rho_coarse_ms.gpl'
   
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    data = madedata[0]
    T = data.size / len(madedata[1]) 		# extent in time dir
    tag01 = madedata[1][0]
    tag02 = madedata[1][1]
    print("time extent is: ", T)
    print(dfile)
    suggestedsvdcut = madedata[2]
    
    pfile = None #"vector_fit.p" # last fit

    tmin = 2
    svdcut = 1e-10 #suggestedsvdcut

    fitter = CorrFitter(models=build_models(tag01,tag02,tmin,T))
    for nexp in [2,3,4,5]:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp),p0=None,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)


    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    tstar = tstr 

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

    print('THIS IS THE INFORMATION CENTRE')
    print(len(newdata[tags[0]]), len(newdata[tags[1]]))
    print('tmin, tstar', tmin, tstar)

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamus = g2.a_mu(vpol,1/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamus = g2.a_mu(vpol,1/3.)


    amus_rt = chargedamus/unchargedamus
    amu_diff = chargedamus-unchargedamus
    amu_diff_rt = amu_diff/unchargedamus

    print('[s] a_mu[QCD+QED]{0:>20}\n|a_mu[QCD+QED]-a_mu[QCD]{1:>20}\n|a_mu/a_mu {2:>20}'.format(chargedamus,amu_diff, amu_diff_rt))
######################################################################################
    
def build_prior(nexp):
    prior = gv.BufferDict()

    prior.add('log(a1:vec:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:s)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:s)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:s)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:s)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:s)'][0] = log(gv.gvar(1,1))
    #prior['logdE:etac'][0] = log(gv.gvar(0.1,0.05))

    prior.add('log(a1:vec:qed:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:qed:s)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:qed:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:qed:s)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:qed:s)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:qed:s)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:qed:s)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:qed:s)'][0] = log(gv.gvar(1,1))

    return prior
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
    return models
##    

def make_data(filename,norm=1):
    dset = gv.dataset.Dataset(filename)
    s = gv.dataset.svd_diagnosis(dset)
    print(dset.keys())
    for tag in dset.keys():
        dset[tag] = norm*np.array(dset[tag])
    return (gv.dataset.avg_data(dset), dset.keys(), s.svdcut)


if __name__ == '__main__':
    for i in range(13,14):
        main(i)
