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
from commonweal import hbarc, w0overa, w0, ZV, ZVqedd, ZVqedu

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'vc'

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

a = (w0/w0overa[a_str])/hbarc

pfile = None #"./fits/mphys15.p" # last fit

store = None #"./store/vcoarse/mphys15.p"

def main(tstr):

    dfile = '/home/gray/Desktop/lattice-analysis/data/qqed/vcoarse/mud_vt_vcoarse.gpl'
    madedata = make_data(dfile,norm=3., binsize=4) # factor of 3 for colour (missed in extraction)
    data = madedata[0]
    T = len(data) / len(madedata[1]) 		# extent in time dir
    T = int(T)
    tag01 = list(madedata[1])[0]
    tag02 = list(madedata[1])[1]
    tag03 = list(madedata[1])[2]		# =tag01 if no isospin breaking
    tag04 = list(madedata[1])[3]
    #sys.exit(0)

    #ratiodata = {}
    #ratiodata['up'] = data[tag02]/data[tag01]
    #ratiodata['down'] = data[tag04]/data[tag03]

    suggestedsvdcut = madedata[2]
    

    tmin = 2
    svdcut = 1e-10 #suggestedsvdcut

    fitter = CorrFitter(models=build_models_phys(tag01,tag02,tag03,tag04, tmin, T))
    for nexp in [2,3,4,5]:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp,3.1,a.mean),p0=pfile,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)


    tdata = range(T)
    #tfit = range(tmin,T+1-tmin) # all ts
    tfit = range(tmin, 15) # all ts
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
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqedu, ainv=1/a, periodic=False)
    chargedamuu = g2.a_mu(vpol,2/3.)

    #moments = g2.moments(newdata[tags[2]],Z=ZV,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamud = g2.a_mu(vpol,1/3.)
 
    #moments = g2.moments(newdata[tags[3]],Z=ZVqed,ainv=1/a,periodic=False)
    #vpol = g2.vacpol(moments,order=(2,1))
    vpol = g2.fourier_vacpol(newdata[tags[3]], Z=ZVqedd, ainv=1/a, periodic=False)
    chargedamud = g2.a_mu(vpol,1/3.)

    amu_qcd = unchargedamud+unchargedamuu
    amu_qcdqed = chargedamuu+chargedamud
    amu_diff = amu_qcdqed-amu_qcd
    d_diff = chargedamud - unchargedamud
    u_diff = chargedamuu - unchargedamuu

    print("down diff is %s"%(d_diff))
    print("up diff is %s"%(u_diff))
    print("amu qcd is %s"%(amu_qcd))
    print("amu diff is %s"%(amu_diff))

    #inputs = { 'Q0':fitdataunchargedup.p['dE:vec:qed:u'], 'a':1/ainv, 'Zv':Zv, 'Zvqed':Zvqed }
    inputs = { 'fit':fit.p['ao:vec:u'], 'a':a, 'Zv':ZV, 'Zvqed':ZVqed }
    outputs = { 'up':unchargedamuu, 'up-diff':u_diff ,'down-qed':chargedamud, 'down-diff':d_diff }
    print(gv.fmt_errorbudget(outputs, inputs, ndecimal=1))

    store_results = {"up-nocharge":unchargedamuu, "down-nocharge":unchargedamud, "up-charge":chargedamuu, "down-charge":chargedamud}

    gv.dump(store_results, store, add_dependencies=True)

############################################################################

if __name__ == '__main__':
    for i in range(9, 10):
        main(i)
