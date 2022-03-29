#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
Edited by Gaurav Ray February 2022
"""

import matplotlib
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
from commonweal import w0overa, w0, ZV, ZVqed_u, ZVqed_d, hbarc
from fitting import make_data, build_prior, build_models_ml

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'f'
w0overa = w0overa[a_str]	

ZV = ZV[a_str] 
ZVqedu = ZVqed_u[a_str]*ZV 
ZVqedd = ZVqed_d[a_str]*ZV 
print("ZVqed\nQ=2/3: %s\nQ=1/3: %s"%(ZVqedu,ZVqedd))
a = (w0/w0overa)/hbarc		# in units of (GeV)^-1

store_fit = './fits/fine/3ml.p'
store_results = None #'./store/coarse/7ml.p'

print("lattice spacing: ", a)

bnsze = 1
tmin = 2
T_STAR = [23]

## tstar reconstruction
store_corr = './fits/tstar_fine_3ml.p'

def main(tstr):
    dfile = '/home/gray/Desktop/lattice-analysis/data/qqed/fine/3ml_rho_finep.gpl'
   
    madedata = make_data(dfile,norm=3.,binsize=bnsze) # factor of 3 for colour (missed in extraction)
    data = gv.dataset.avg_data(madedata[0])
    tag01 = list(madedata[1])[0]
    tag02 = list(madedata[1])[1]
    tag03 = list(madedata[1])[2]	
    T = data.size / len(madedata[1]) 		# extent in time dir
    T = int(T)
    print("T = %d\ntfit=[%d,%d]"%(T,tmin,T-tmin))
   ### svd diagnosis ###
    s = gv.dataset.svd_diagnosis(madedata[0], models=build_models_ml(tag01,tag02,tag03,tmin,T))
    #s.plot_ratio(show=True)
    ####################
    suggestedsvdcut = s.svdcut
    svdcut = suggestedsvdcut
    print('svd cut = %f'%svdcut)

    fitter = CorrFitter(models=build_models_ml(tag01,tag02,tag03,tmin,T))
    prior = build_prior(4,5,a.mean)

    p0 = store_fit
    for nexp in [2,3,4]:
        print("N=%d"%nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0,maxit=50000,svdcut=svdcut,add_svdnoise=False, nterm=nexp)
        p0=fit.pmean
    print(fit.format(pstyle='m'))

#    ### NOISY FIT ###
#    noisyfit = fitter.lsqfit(data=data,prior=build_prior(5,2.0,a.mean),p0=store_fit,maxit=20000,svdcut=svdcut,add_svdnoise=True)
#    print("ADD SVD NOISE FIT")
#    print(noisyfit.format(pstyle='m'))
#    ##################
    
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    tstar = tstr 

    tcut = 200

    newdata = {}

    tags = ['nocharge', 'up-charge', 'down-charge']
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

    print(len(newdata[tags[0]]), len(newdata[tags[1]]), len(newdata[tags[2]]))
    print('tmin, tstar', tmin, tstar)

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
    unchargedamuu = g2.a_mu(vpol,2/3.)
    unchargedamud = g2.a_mu(vpol,1/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqedu, ainv=1/a, periodic=False)
    chargedamuu = g2.a_mu(vpol,2/3.)
 
    vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZVqedd, ainv=1/a, periodic=False)
    chargedamud = g2.a_mu(vpol,1/3.)

    amu_qcd = unchargedamud+unchargedamuu
    amu_qcdqed = chargedamuu+chargedamud
    amu_rt = amu_qcdqed/amu_qcd
    amu_diff = amu_qcdqed-amu_qcd

    print("amu WITHOUT qed is %s"%(amu_qcd))
    print("amu WITH qed is %s"%(amu_qcdqed))
    print("QED+QCD/QCD is %s"%(amu_rt))
    print("[QED+QCD]-QCD = %s"%(amu_diff))

    to_save = {"up-nocharge":unchargedamuu, "up-charge":chargedamuu, "down-nocharge":unchargedamud, "down-charge":chargedamud, "mq":[ amu_qcd, amu_diff, amu_rt ] }
    corr = {"corr":data, "newcorr":newdata}

    gv.dump( to_save, store_results, add_dependencies=True ) 
    gv.dump( corr, store_corr, add_dependencies=True )
######################################################################################
    
if __name__ == '__main__':
    for i in T_STAR:
        main(i)
