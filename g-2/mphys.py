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
from fitting import build_prior, make_data, build_models_phys
sys.path.append('./plotters')
from commonweal import hbarc, w0overa, w0, ZV, ZVqed

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'vc'

ZV = ZV[a_str]
ZVqed = ZVqed[a_str]*ZV

a = (w0/w0overa[a_str])/hbarc

pfile = "./stored-fits/vector_mud.p" # last fit

def main(tstr):

    dfile = '/home/gray/Desktop/lattice-analysis/data/qqed/vcoarse/mphys_vcoarse_full.gpl'
#    tag01 = 'rho_m' + str(mq)
#    tag02 = 'rho_m' + str(mq) + '_ucav'
#    tag03 = 'rho_m' + str(mq)
#    tag04 = 'rho_m' + str(mq) + '_dcav'
   
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    data = madedata[0]
    T = data.size / len(madedata[1]) 		# extent in time dir
    tag01 = madedata[1][0]
    tag02 = madedata[1][1]
    tag03 = madedata[1][2]		# =tag01 if no isospin breaking
    tag04 = madedata[1][3]
    #sys.exit(0)

    #ratiodata = {}
    #ratiodata['up'] = data[tag02]/data[tag01]
    #ratiodata['down'] = data[tag04]/data[tag03]

    suggestedsvdcut = madedata[2]
    

    tmin = 2
    svdcut = 1e-10 #suggestedsvdcut

    fitter = CorrFitter(models=build_models_phys(tag01,tag02,tag03,tag04, tmin, T))
    for nexp in [2,3,4,5]:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp,1.5,a.mean),p0=pfile,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)


    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    tstar = tstr 

    tcut = 200

    newdata = {}

    tags = ['up', 'up_qed', 'down', 'down_qed']
    newdata[tags[0]] = data[tag01][:tstar+1]
    newdata[tags[1]] = data[tag02][:tstar+1]
    newdata[tags[2]] = data[tag03][:tstar+1]
    newdata[tags[3]] = data[tag04][:tstar+1]


    # do replacement of data with fit

    fitdataunchargedup = Corr2(datatag=tags[0],tdata=range(T),a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)).fitfcn(fit.p)
    fitdataunchargeddown = Corr2(datatag=tags[2],tdata=range(T),a=('a1:vec:d','ao:vec:d'),b=('a1:vec:d','ao:vec:d'),dE=('dE:vec:d','dEo:vec:d'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(48):
            if index >= tcut:
                newdata[tags[0]] = np.append(newdata[tags[0]],0.)
                newdata[tags[2]] = np.append(newdata[tags[2]],0.)
            elif index > tstar:
                newdata[tags[0]] = np.append(newdata[tags[0]],fitdataunchargedup[index])
                newdata[tags[2]] = np.append(newdata[tags[2]],fitdataunchargeddown[index])

    fitdatachargedup = Corr2(datatag=tags[1],tdata=range(T),a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)).fitfcn(fit.p)
    fitdatachargedown = Corr2(datatag=tags[3],tdata=range(T),a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(48):
            if index >= tcut:
                newdata[tags[1]] = np.append(newdata[tags[1]],0.)
                newdata[tags[3]] = np.append(newdata[tags[3]],0.)
            elif index > tstar:
                newdata[tags[1]] = np.append(newdata[tags[1]],fitdatachargedup[index])
                newdata[tags[3]] = np.append(newdata[tags[3]],fitdatachargedown[index])

    print('THIS IS THE INFORMATION CENTRE')
    print(len(newdata[tags[0]]), len(newdata[tags[1]]), len(newdata[tags[2]]), len(newdata[tags[3]]))
    print('tmin, tstar', tmin, tstar)

    #moments = g2.moments(newdata[tags[0]],Z=ZV,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamuu = g2.a_mu(vpol, 2/3.)
 
    #moments = g2.moments(newdata[tags[1]],Z=ZVqed,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamuu = g2.a_mu(vpol,2/3.)

    #moments = g2.moments(newdata[tags[2]],Z=ZV,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamud = g2.a_mu(vpol,1/3.)
 
    #moments = g2.moments(newdata[tags[3]],Z=ZVqed,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[3]], Z=ZVqed, ainv=1/a, periodic=False)
    chargedamud = g2.a_mu(vpol,1/3.)

    d_rt = chargedamud/unchargedamud
    u_rt = chargedamuu/unchargedamuu
    amu_qcd = unchargedamud+unchargedamuu
    amu_qcdqed = chargedamuu+chargedamud
    amu_rt = amu_qcdqed/amu_qcd
    amu_diff = amu_qcdqed-amu_qcd

    d_diff = chargedamud - unchargedamud
    u_diff = chargedamuu - unchargedamuu

    print("down diff is %s"%(d_diff))
    print("up diff is %s"%(u_diff))

    store_results = {"qcd":amu_qcd, "qcd+qed":amu_qcdqed, "rt":amu_rt, "diff":amu_diff}

    gv.dump(store_results, './results/mphys.p')

############################################################################

if __name__ == '__main__':
    for i in range(13, 14):
        main(i)
