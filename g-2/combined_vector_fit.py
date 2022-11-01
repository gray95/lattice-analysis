#!/usr/bin/env python
# encoding: utf-8
"""

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
Edited by Dan Hatton December 2016
Edited by Gaurav Ray April 2022
"""

import os
import sys
import dill as pickle
import bz2
#---------------------#
import numpy as np
import corrfitter as cf
from corrfitter import Corr2, CorrFitter
import gvar as gv
import g2tools as g2
import lsqfit as lsq
import collections
#---------------------#
sys.path.append('./plotters')
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
#from commonweal import base, fname_vt
from fitting import make_data

############################################################################

lsq.LSQFit.fmt_parameter = '%8.4f +- %8.4f'
a_str = 'f' 

a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))
ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

#name = fname_vt[a_str]
#datafile = os.path.join(base, name)
datafile="/home/gray/Desktop/lattice-analysis/data/qqed/ms-tuning-fine/mistune_rho_fine.gpl"
masses = ['ms-1', 'ms', 'ms1', 'ms2']
bnSze = 1
tmin = 2
NEXP = [1,2,3,4,5]
Nmax = 7
# set t* ~ 2fm
T_STAR = [int(np.ceil(2.5/(a.mean*hbarc)))]

store_fitp = './fits/post/vt_mistuned_'+a_str+'_bin'+str(bnSze)+'.p'
store_res  = None #'./fits/vt_'+a_str+'_bin'+str(bnSze)+'.p'
store_FIT = None #'./fits/fitobj/'+a_str+'_FIT_bin'+str(bnSze)+'.pbz2'

############################################################################


def main(tstr):
    dfile = datafile
   
    madedata = make_data(dfile,norm=3.,binsize=bnSze) # factor of 3 for colour (missed in extraction)
    data = gv.dataset.avg_data(madedata[0])
    tags = data.keys()
    T = data.size / len(tags) 		# extent in time dir
    T = int(T)
    print("time extent = %d"%(T))
    print("t* = %d = %sfm"%(tstr, a*hbarc*tstr))
    TAGS = {'Q0':[s for s in tags if 'q0' in s], 'Q1':[s for s in tags if 'q1' in s], 'Q2':[s for s in tags if 'q2' in s]}
    ### svd diagnosis ###
    s = gv.dataset.svd_diagnosis(madedata[0], models=make_models(masses, TAGS, T, tmin))
#    s.plot_ratio(show=True)
    svdcut = s.svdcut
    ####################

    print('svd cut = %f'%svdcut)
    fitter = CorrFitter(models=make_models(masses, TAGS, T, tmin))
    NMASSES = len(masses)
    p0 = store_fitp	# saved fit priors 
    prior = make_prior(Nmax,masses,a.mean)
    for nexp in NEXP:
        print("nterm = %d"%nexp)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, nterm=nexp, maxit=30000,svdcut=svdcut, debug=True)
        p0 = fit.pmean

    print(fit.format(pstyle='m'))

    ### NOISY FIT ###
    #noisyfit = fitter.lsqfit(data=data,prior=prior, p0=p0, nterm=NEXP[-1], maxit=30000, svdcut=svdcut, noise=(True,False))
    #print("ADDING SVD NOISE TO FIT\n######################")
    #print(noisyfit.format(pstyle='m'))
    ##################

    tcut = T

    tstar = tstr 

    res = gv.BufferDict()
    res["down-nocharge"] = gv.gvar(np.zeros(4))
    res["down-charge"]   = gv.gvar(np.zeros(4))
    res["up-nocharge"]   = gv.gvar(np.zeros(4))
    res["up-charge"]     = gv.gvar(np.zeros(4))

    newdata = gv.BufferDict()

    # do replacement of data with fit
    for M in range(NMASSES):
        m=masses[M]+':'
        tags = ['no-charge', 'down-charge', 'up-charge']

        fitdatauncharged = Corr2(datatag=tags[0],tdata=range(tcut),a=(m+'a:n',m+'ao:n'),b=(m+'a:n',m+'ao:n'),dE=(m+'dE:n',m+'dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
        fitdatachargeddown = Corr2(datatag=tags[1],tdata=range(tcut),a=(m+'a:d',m+'ao:d'),b=(m+'a:d',m+'ao:d'),dE=(m+'dE:d',m+'dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
        fitdatachargedup = Corr2(datatag=tags[2],tdata=range(tcut),a=(m+'a:u',m+'ao:u'),b=(m+'a:u',m+'ao:u'),dE=(m+'dE:u',m+'dEo:u'),s=(1.,-1.)).fitfcn(fit.p)

        newdata[tags[0]] = np.append( data[TAGS['Q0'][M]][:tstar+1], fitdatauncharged[tstar+1:] )
        newdata[tags[1]] = np.append( data[TAGS['Q1'][M]][:tstar+1], fitdatachargeddown[tstar+1:] )
        newdata[tags[2]] = np.append( data[TAGS['Q2'][M]][:tstar+1], fitdatachargedup[tstar+1:] )

        # fourier transform to find vacuum pol; intergrate vpol over kernel for amu
        vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, ainv=1/a, periodic=False)
        unchargedamud = g2.a_mu(vpol,1/3.)
        unchargedamuu = g2.a_mu(vpol,2/3.)
     
        vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqedd, ainv=1/a, periodic=False)
        chargedamud = g2.a_mu(vpol,1/3.)

        vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZVqedu, ainv=1/a, periodic=False)
        chargedamuu = g2.a_mu(vpol,2/3.)

        print("MASS #%d"%M)
        print("strange quark")
        amu = unchargedamud
        amuqed = chargedamud

        amu_diff = amuqed-amu
        amu_rt = amuqed/amu
        print("amu is %s"%amu)
        print("amu with qed is %s"%amuqed)
        print("diff is %s"%amu_diff)

        res['down-nocharge'][M] = unchargedamud
        res['up-nocharge'][M]   = unchargedamuu
        res['down-charge'][M]   = chargedamud
        res['up-charge'][M]     = chargedamuu

        # error analysis
        modeldata = gv.BufferDict()
        corrdata  = gv.BufferDict()
        for k in newdata:
            modeldata[k] = newdata[k][tstar+1:]
            corrdata[k]  = newdata[k][:tstar+1]

        inputs = { "lat_spacing":a, "ZV":ZV, "data":fit.y, "model":modeldata,
		   "model+data":newdata, "svdcut":fit.correction }
        outputs = { 'amu':amu, 'diff':amu_diff, 'rt':amu_rt }
        print(gv.fmt_errorbudget(outputs=outputs, inputs=inputs))

#    gv.dump( res, store_res,  add_dependencies=True)

#    with bz2.BZ2File(store_FIT, 'wb') as f:
#        pickle.dump( fit, f )


    ######################################################################################

def make_models(masses, keys, T, tmin, otherkeys=None):
    """ Create corrfitter model for G(t). """
    TDATA = range(T)
    TP = T
    TFIT = range(tmin,T+1-tmin)
    models = []
    for (k0,k1,k2,m) in zip(keys['Q0'],keys['Q1'],keys['Q2'],masses):
        c = m+':'
        models.append(cf.Corr2(datatag=k0, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:n',c+'ao:n'), b=(c+'a:n',c+'ao:n'), dE=(c+'dE:n',c+'dEo:n'), s=(1.,-1.)))
        models.append(cf.Corr2(datatag=k1, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:d',c+'ao:d'), b=(c+'a:d',c+'ao:d'), dE=(c+'dE:d',c+'dEo:d'), s=(1.,-1.)))
        models.append(cf.Corr2(datatag=k2, otherdata=otherkeys, tdata=TDATA, tfit=TFIT, tp=TP, a=(c+'a:u',c+'ao:u'), b=(c+'a:u',c+'ao:u'), dE=(c+'dE:u',c+'dEo:u'), s=(1.,-1.)))
    return models

def make_prior(N, MASSES, a):
    """ Create prior for N-state fit. """
    prior = gv.BufferDict()
    for m in MASSES:
        prior['log('+m+':a:n)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':a:u)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':a:d)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':dE:n)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]
        prior['log('+m+':dE:u)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]
        prior['log('+m+':dE:d)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]

        prior['log('+m+':ao:n)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':ao:u)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':ao:d)'] = gv.log(gv.gvar(['0.1(1)'] + (N-1)*['0.01(99)']))
        prior['log('+m+':dEo:n)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]
        prior['log('+m+':dEo:u)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]
        prior['log('+m+':dEo:d)'] = [gv.log(gv.gvar(a,a)) for i in range(N)]
    return prior

def print_results(fit, MASSES):
    p = fit.p
    #g = collections.OrderedDict()
    g = gv.BufferDict()
    for i in MASSES:
        m = i+':'
        for Q in ['n', 'u', 'd']:
            g[m+'E:'+Q] = np.cumsum(p[m+'dE:'+Q])[0]
            g[m+'a:'+Q] = p[m+'a:'+Q][0]
    print(gv.tabulate(g, ncol=2))
###############################################################################    

if __name__ == '__main__':
    for i in T_STAR:
        main(i)
