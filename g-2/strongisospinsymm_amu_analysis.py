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
sys.path.append('./plotters')
from corrfitter import Corr2,Corr3,CorrFitter
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import math
from math import exp as num_exp
import matplotlib.pyplot as plt
import g2tools as g2
from commonweal import w0overa, w0, ZV, ZVqed

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'c'
w0overa = w0overa[a_str]	
ZV = ZV[a_str] 
ZVqed = ZVqed[a_str]*ZV 

hbarc = 0.197326968
a = (w0/w0overa)/hbarc		# in units of (GeV)^-1

print("lattice spacing: ", (w0/w0overa))

def main(tstr):
    dfile = '/home/gray/Desktop/lattice-analysis/data/qqed/coarse/3ml_rho_coarse.gpl'
   
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    data = madedata[0]
    T = data.size / len(madedata[1]) 		# extent in time dir
    tag01 = madedata[1][0]
    tag02 = madedata[1][1]
    tag03 = madedata[1][2]		# =tag01 if no isospin breaking
    print("time extent is: ", T)
    print(dfile)
    suggestedsvdcut = madedata[2]
    
    pfile = None #"vector_fit.p" # last fit

    tmin = 2
    svdcut = suggestedsvdcut

    fitter = CorrFitter(models=build_models(tag01,tag02,tag03,tmin,T))
    for nexp in [2,3,4,5]:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp),p0=None,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)


    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    tstar = tstr 

    tcut = 200

    newdata = {}

    tags = ['nocharge', 'up_charge', 'down-charge']
    newdata[tags[0]] = data[tag01][:tstar+1]
    newdata[tags[1]] = data[tag02][:tstar+1]
    newdata[tags[2]] = data[tag03][:tstar+1]


    # do replacement of data with fit

    fitdatauncharged = Corr2(datatag=tags[0],tdata=range(T),a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(T):
            if index >= tcut:
                newdata[tags[0]] = np.append(newdata[tags[0]],0.)
            elif index > tstar:
                newdata[tags[0]] = np.append(newdata[tags[0]],fitdatauncharged[index])

    fitdatachargedup = Corr2(datatag=tags[1],tdata=range(T),a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)).fitfcn(fit.p)
    fitdatachargedown = Corr2(datatag=tags[2],tdata=range(T),a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(T):
            if index >= tcut:
                newdata[tags[1]] = np.append(newdata[tags[1]],0.)
                newdata[tags[2]] = np.append(newdata[tags[2]],0.)
            elif index > tstar:
                newdata[tags[1]] = np.append(newdata[tags[1]],fitdatachargedup[index])
                newdata[tags[2]] = np.append(newdata[tags[2]],fitdatachargedown[index])

    print('THIS IS THE INFORMATION CENTRE')
    print(len(newdata[tags[0]]), len(newdata[tags[1]]), len(newdata[tags[2]]))
    print('tmin, tstar', tmin, tstar)

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamuu = g2.a_mu(vpol,2/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamuu = g2.a_mu(vpol,2/3.)

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamud = g2.a_mu(vpol,1/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamud = g2.a_mu(vpol,1/3.)

    d_rt = chargedamud/unchargedamud
    u_rt = chargedamuu/unchargedamuu 
    amu_qcd = unchargedamud+unchargedamuu
    amu_qcdqed = chargedamuu+chargedamud
    amu_rt = amu_qcdqed/amu_qcd
    amu_diff_rt = (amu_qcdqed-amu_qcd)/amu_qcd
    amu_diff = amu_qcdqed-amu_qcd

    print('[d] a_mu[QCD+QED]||a_mu[QCD+QED]/a_mu[QCD]{0:>20}||{1}'.format(chargedamud,d_rt))
    print('[u] a_mu[QCD+QED]||a_mu[QCD+QED]/a_mu[QCD]{0:>20}||{1}'.format(chargedamuu,u_rt))
    print('\n[u+d] a_mu[QCD+QED]||ratio||diff||diff rt\n{0:20}{1:20}{2:20}{3:20}'.format(amu_qcdqed, amu_rt, amu_diff, amu_diff_rt))
######################################################################################
    
def build_prior(nexp):
    prior = gv.BufferDict()

    prior.add('log(a1:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:u)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(ao:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:u)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(dE:vec:u)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:u)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:u)'][0] = log(gv.gvar(1,1))

    prior.add('log(a1:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:qed:u)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(ao:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:qed:u)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(dE:vec:qed:u)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:qed:u)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:qed:u)'][0] = log(gv.gvar(1,1))

    prior.add('log(a1:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:d)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(ao:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:d)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(dE:vec:d)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:d)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:d)'][0] = log(gv.gvar(1,1))

    prior.add('log(a1:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:qed:d)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(ao:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:qed:d)'][0] = log(gv.gvar(0.5,0.5))
    prior.add('log(dE:vec:qed:d)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:qed:d)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:qed:d)'][0] = log(gv.gvar(1,1))

    return prior
##

def build_models(tag01,tag02,tag03,tmin,T):
 
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    
    tp = T # periodic
    models = [
        Corr2(datatag=tag01,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),
            dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)),
        Corr2(datatag=tag02,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),
            dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)),
        Corr2(datatag=tag03,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),
            dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.))
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
    for i in range(17,18):
        main(i)