#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
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

DISPLAYPLOTS = False
#tmin = 8
#svdcut = 1.0e-4
#svdcut = 0.003

hbarc = 0.197326968

def main():
    dfile = 'rho_vcphys_bothcharges_up_down.gpl'
    tag01 = 'uncharged-up'
    tag02 = 'charged-up'
    tag03 = 'uncharged-down'
    tag04 = 'charged-down'

    T = 48
    madedata = make_data(dfile,norm=3.) # factor of 3 for colour (missed in extraction)
    data = madedata[0]

    ratiodata = {}
    ratiodata['up'] = data['charged-up']/data['uncharged-up']
    ratiodata['down'] = data['charged-down']/data['uncharged-down']

    suggestedsvdcut = madedata[1]
    
    pfile = "vector_fit.p" # last fit

    tmin = 2
    svdcut = 1e-16


    fitter = CorrFitter(models=build_models(tag01,tag02,tag03,tag04,tmin,T))
    for nexp in [2,3,4,5]:
        fit = fitter.lsqfit(data=data,prior=build_prior(nexp),p0=pfile,maxit=20000,svdcut=svdcut,add_svdnoise=False)
    print(fit)


    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T

    tstar = 14

    tcut = 200

    newdata = {}

    newdata['uncharged-up'] = data['uncharged-up'][:tstar]
    newdata['charged-up'] = data['charged-up'][:tstar]
    newdata['uncharged-down'] = data['uncharged-down'][:tstar]
    newdata['charged-down'] = data['charged-down'][:tstar]


    # do replacement of data with fit

    fitdataunchargedup = Corr2(datatag='uncharged-up',tdata=range(48),a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)).fitfcn(fit.p)
    fitdataunchargeddown = Corr2(datatag='uncharged-down',tdata=range(48),a=('a1:vec:d','ao:vec:d'),b=('a1:vec:d','ao:vec:d'),dE=('dE:vec:d','dEo:vec:d'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(48):
            if index >= tcut:
                newdata['uncharged-up'] = np.append(newdata['uncharged-up'],0.)
                newdata['uncharged-down'] = np.append(newdata['uncharged-down'],0.)
            elif index >= tstar:
                newdata['uncharged-up'] = np.append(newdata['uncharged-up'],fitdataunchargedup[index])
                newdata['uncharged-down'] = np.append(newdata['uncharged-down'],fitdataunchargeddown[index])

    fitdatachargedup = Corr2(datatag='charged-up',tdata=range(48),a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)).fitfcn(fit.p)
    fitdatachargedown = Corr2(datatag='charged-down',tdata=range(48),a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.)).fitfcn(fit.p)
    for index in range(48):
            if index >= tcut:
                newdata['charged-up'] = np.append(newdata['charged-up'],0.)
                newdata['charged-down'] = np.append(newdata['charged-down'],0.)
            elif index >= tstar:
                newdata['charged-up'] = np.append(newdata['charged-up'],fitdatachargedup[index])
                newdata['charged-down'] = np.append(newdata['charged-down'],fitdatachargedown[index])

    w0 = gv.gvar('0.1715(9)')
    w0overa = gv.gvar('1.1367(5)')
    ##ZV = gv.gvar('0.9881(10)')
    ZV = gv.gvar('0.9837(20)')
    ZVqed = gv.gvar('0.999544(14)')*ZV
    hbarc = 0.197326968
    a = (w0/w0overa)/hbarc

    moments = g2.moments(newdata['uncharged-up'],Z=ZV,ainv=1/a,periodic=False)
    vpol = g2.vacpol(moments,order=(2,1))
    unchargedamuu = g2.a_mu(vpol,2/3.)

    print('up a_mu[QCD]: ',unchargedamuu)
 
    moments = g2.moments(newdata['charged-up'],Z=ZVqed,ainv=1/a,periodic=False)
    vpol = g2.vacpol(moments,order=(2,1))
    chargedamuu = g2.a_mu(vpol,2/3.)

    print('up a_mu[QCD+QED]: ',chargedamuu)

    moments = g2.moments(newdata['uncharged-down'],Z=ZV,ainv=1/a,periodic=False)
    vpol = g2.vacpol(moments,order=(2,1))
    unchargedamud = g2.a_mu(vpol,1/3.)

    print('down a_mu[QCD]: ',unchargedamud)
 
    moments = g2.moments(newdata['charged-down'],Z=ZVqed,ainv=1/a,periodic=False)
    vpol = g2.vacpol(moments,order=(2,1))
    chargedamud = g2.a_mu(vpol,1/3.)

    print('down a_mu[QCD+QED]: ',chargedamud)

    print('up a_mu[QCD+QED]/a_mu[QCD]: ',chargedamuu/unchargedamuu)
    #print(chargedamud/unchargedamud)
    print('up+down a_mu[QCD]: ',unchargedamud+unchargedamuu)
    print('up+down a_mu[QCD+QED]/a_mu[QCD]: ',(chargedamud+chargedamuu)/(unchargedamud+unchargedamuu))



    
##

    
def build_prior(nexp):
    prior = gv.BufferDict()

    prior.add('log(a1:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:u)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:u)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:u)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:u)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:u)'][0] = log(gv.gvar(1,15))
    #prior['logdE:etac'][0] = log(gv.gvar(0.1,0.05))

    prior.add('log(a1:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:qed:u)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:qed:u)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dE:vec:qed:u)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:qed:u)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:qed:u)'][0] = log(gv.gvar(1,1))

    prior.add('log(a1:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:d)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:d)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:d)',[log(gv.gvar(2,2)) for i in range(nexp)])
    prior['log(dE:vec:d)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:d)'][0] = log(gv.gvar(1,15))
    #prior['logdE:etac'][0] = log(gv.gvar(0.1,0.05))

    prior.add('log(a1:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(a1:vec:qed:d)'][0] = log(gv.gvar(0.5,1))
    prior.add('log(ao:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(ao:vec:qed:d)'][0] = log(gv.gvar(0.5,1))
    #prior.add('as1:etac',[gv.gvar(0.001,0.01) for i in range(nexp)])
    prior.add('log(dE:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dE:vec:qed:d)'][0] = log(gv.gvar(1,1))
    prior.add('log(dEo:vec:qed:d)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dEo:vec:qed:d)'][0] = log(gv.gvar(1,1))

    return prior
##

def build_models(tag01,tag02,tag03,tag04,tmin,T):
 
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
            a=('a1:vec:d','ao:vec:d'),b=('a1:vec:d','ao:vec:d'),
            dE=('dE:vec:d','dEo:vec:d'),s=(1.,-1.)),
        Corr2(datatag=tag04,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),
            dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.))
    ]
    return models
##    

def make_data(filename,norm=1):
    dset = gv.dataset.Dataset(filename)
    s = gv.dataset.svd_diagnosis(dset)
    print(s.svdcut)
    for tag in dset.keys():
        dset[tag] = norm*np.array(dset[tag])
    #print(len(dset['etac:l.l']))
    #print len(dset['etac_ptpt_smsm_pair00'])
    return (gv.dataset.avg_data(dset),s.svdcut)


if __name__ == '__main__':
    main()

