#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
Edited by Gaurav Ray March 2021
"""

import os
import sys
sys.path.append('./plotters')
import corrfitter as cf
from corrfitter import Corr2,Corr3,CorrFitter
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import math as m
from math import exp as num_exp
import g2tools as g2
import lsqfit as lsq
import collections
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
from fitting import make_data

lsqfit.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'vc' 

a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

base = 'vcoarse/'
name = 'Nml_rho_vc'
datafile = os.path.join('../data/qqed', base, name+'.gpl')
masses = ['3ml', '5ml', '7ml']
bnSze = 4
tmin = 2
T_STAR = [13]

store_fit = './fits/vt_vc.p'
store_res = './store/vc.p'
store_corr= './fits/vc_rhocorr.p'

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

    TAGS = {'Q0':['VTvc3mlq0','VTvc5mlq0','VTvc7mlq0'], 'Q101':['VTvc3mlq1','VTvc5mlq1','VTvc7mlq1'], 'Q202':['VTvc3mlq2','VTvc5mlq2','VTvc7mlq2']}

    ### svd diagnosis ###
    s = gv.dataset.svd_diagnosis(madedata[0], models=make_models(TAGS, T, tmin))
#    s.plot_ratio(show=True)
    suggestedsvdcut = s.svdcut
    svdcut = suggestedsvdcut
    ####################

    print('svd cut = %f'%svdcut)
    fitter = CorrFitter(models=make_models(keys=TAGS, T=T, tmin=tmin))
    NEXP = [1,2,3,4,5]
    NMASSES = len(masses)
    p0 = store_fit
    for nexp in NEXP:
        print("nterm = %d"%nexp)
        fit = fitter.lsqfit(data=data,prior=make_prior(nexp, NMASSES), p0=p0 ,maxit=20000,svdcut=svdcut,add_svdnoise=False)
        p0 = fit.pmean

    tstar = tstr 
    print(fit.format(pstyle='m'))

    #print_results(fit, ainv, ZV, ZVqed)

    tcut = 200

    newdata = {}

    res = {"diff":[], "ratio":[], "amu":[], "amuqed":[], 
            'down-nocharge':[], 'up-nocharge':[], 'down-charge':[], 'up-charge':[]}
    corr= {}

    # do replacement of data with fit
    for M in range(NMASSES):
        m='m'+str(M)
        tags = ['nocharge', 'down-charge', 'up-charge']
        newdata[tags[0]] = data[TAGS['Q0'][M]][:tstar+1]
        newdata[tags[1]] = data[TAGS['Q101'][M]][:tstar+1]
        newdata[tags[2]] = data[TAGS['Q202'][M]][:tstar+1]


        fitdatauncharged = Corr2(datatag=tags[0],tdata=range(T),a=(m+'a:n',m+'ao:n'),b=(m+'a:n',m+'ao:n'),dE=(m+'dE:n',m+'dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
        for index in range(T):
                if index >= tcut:
                    newdata[tags[0]] = np.append(newdata[tags[0]],0.)
                elif index > tstar:
                    newdata[tags[0]] = np.append(newdata[tags[0]],fitdatauncharged[index])

        fitdatachargeddown = Corr2(datatag=tags[1],tdata=range(T),a=(m+'a:d',m+'ao:d'),b=(m+'a:d',m+'ao:d'),dE=(m+'dE:d',m+'dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
        for index in range(T):
                if index >= tcut:
                    newdata[tags[1]] = np.append(newdata[tags[1]],0.)
                elif index > tstar:
                    newdata[tags[1]] = np.append(newdata[tags[1]],fitdatachargeddown[index])

        fitdatachargedup = Corr2(datatag=tags[2],tdata=range(T),a=(m+'a:u',m+'ao:u'),b=(m+'a:u',m+'ao:u'),dE=(m+'dE:u',m+'dEo:u'),s=(1.,-1.)).fitfcn(fit.p)
        for index in range(T):
                if index >= tcut:
                    newdata[tags[2]] = np.append(newdata[tags[2]],0.)
                elif index > tstar:
                    newdata[tags[2]] = np.append(newdata[tags[2]],fitdatachargedup[index])
        corr[masses[M]] = newdata.copy() 
        # fourier transform to find vacuum pol; intergrate vpol over kernel for amu

        vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
        unchargedamud = g2.a_mu(vpol,1/3.)
        unchargedamuu = g2.a_mu(vpol,2/3.)
     
        vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqedd, ainv=1/a, periodic=False)
        chargedamud = g2.a_mu(vpol,1/3.)

        vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZVqedu, ainv=1/a, periodic=False)
        chargedamuu = g2.a_mu(vpol,2/3.)

        print("MASS #%d"%M)
        amu = unchargedamud + unchargedamuu
        amuqed = chargedamud + chargedamuu
        amu_diff = amuqed-amu
        amu_rt = amuqed/amu
        print("amu is %s"%amu)
        print("amu with qed is %s"%amuqed)
        print("diff is %s"%amu_diff)

        res['diff'].append(amu_diff)
        res['ratio'].append(amu_rt)
        res['amu'].append(amu)
        res['amuqed'].append(amuqed)
        res['down-nocharge'].append(unchargedamud)
        res['up-nocharge'].append(unchargedamuu)
        res['down-charge'].append(chargedamud)
        res['up-charge'].append(chargedamuu)

        # error analysis
        inputs = { 'fit':fit.p, 'a':a, 'ZV':ZV, 'ZVqedd':ZVqedd, 'ZVqedu':ZVqedu }
        outputs = { 'amu':amu, 'diff':amu_diff, 'rt':amu_rt }
        print(gv.fmt_errorbudget(outputs=outputs, inputs=inputs))

    corr['org']=data
    gv.dump( res, store_res, True)
    gv.dump( corr, store_corr, True)


    ######################################################################################

def make_models(keys, T, tmin, otherkeys=None):
    """ Create corrfitter model for G(t). """
    TDATA = range(T)
    TP = T
    TFIT = range(tmin,T+1-tmin)
    models = []
    count = 0
    for k in keys['Q0']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:n',c+'ao:n'), b=(c+'a:n',c+'ao:n'), dE=(c+'dE:n',c+'dEo:n'), s=(1.,-1.)))
        count = count + 1
    count = 0
    for k in keys['Q101']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:d',c+'ao:d'), b=(c+'a:d',c+'ao:d'), dE=(c+'dE:d',c+'dEo:d'), s=(1.,-1.)))
        count = count + 1
    count = 0
    for k in keys['Q202']:
        c = 'm'+str(count)
        models.append(cf.Corr2(datatag=k, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:u',c+'ao:u'), b=(c+'a:u',c+'ao:u'), dE=(c+'dE:u',c+'dEo:u'), s=(1.,-1.)))
        count = count + 1
    return models

def make_prior(N, NMASSES):
    """ Create prior for N-state fit. """
    prior = collections.OrderedDict()
    for m in range(NMASSES):
        m = 'm'+str(m)
        prior['log('+m+'a:n)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'a:u)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'a:d)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'dE:n)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dE:u)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dE:d)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))

        prior['log('+m+'ao:n)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'ao:u)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'ao:d)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(0.99)']))
        prior['log('+m+'dEo:n)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dEo:u)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
        prior['log('+m+'dEo:d)'] = gv.log(gv.gvar(['1.0(0.9)'] + (N-1)*['0.5(5)']))
    return prior

def print_results(fit, NMASSES):
    p = fit.p
    g = collections.OrderedDict()
    for i in range(NMASSES):
        m = 'm'+str(i)
        for Q in ['n', 'u', 'd']:
            g[m+'E:'+Q] = np.cumsum(p[m+'dE:'+Q])[0]
            g[m+'a:'+Q] = p[m+'a:'+Q][0]
    print(gv.tabulate(g, ncol=2))
###############################################################################    

if __name__ == '__main__':
    for i in T_STAR:
        main(i)
