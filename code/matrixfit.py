from __future__ import print_function   # makes this work for python2 and 3

import os
import sys
import dill as pickle
import bz2
#---------------------#
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import datetime
from parameters import corrpath, l, key, tmin, tmax, D
from parameters import s_coeff, t0, c_hack, NEXP, bin_size, e_str
from parameters import OSC, NOISE, WRITE_LOG, EIG_SRCs, EIG
from parameters import SRCs, KEYFMT, tag, otag, SAVE_RES, save_res_name, save_fit
sys.path.append('./plotters')
from commonweal import ainv
print(key)
data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=key))     
T = data[key].size
tp=T
TDATA = range(T)
TFIT  = list(TDATA[tmin:tmax+1])
#TFIT.remove(5)
print('t0: %d'%t0, file=l)
print('Temporal extent: %d'%T, file=l)
print('Bin Size: %d'%bin_size, file=l)
print('Fit range %s'%str(TFIT), file=l)
print('DIAG Fit range %s'%str(TFIT[:D]), file=l)

def main():
    print("Fitting:[%s]"%os.path.split(corrpath)[-1], file=l)
    data = make_data(corrpath)
    ## svd cut
    s = gv.dataset.svd_diagnosis(cf.read_dataset(corrpath,binsize=bin_size), models=make_models())
    SVDCUT = s.svdcut
    print("svd cut is: ", SVDCUT, file=l)
    basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRCs, t=(t0, t0+2), tdata=TDATA, osc=OSC)
    if EIG:
        #data = basis.apply(data, keyfmt=(tag+'{s1}',otag+'{s2}'))
        data = basis.apply(data, keyfmt=KEYFMT)
        fitter = cf.CorrFitter(models=make_models_eig())
    else:
        fitter = cf.CorrFitter(models=make_models())

    p0 = None
    for N in NEXP:
        print(30 * '=', 'nterm =', N, file=l)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        if WRITE_LOG:
            print(fit.format(pstyle='m'),file=l)
        p0 = fit.pmean
    
    if fit.Q < 0.1:
        print(" BAD FIT - ABORTING ")
        print("Q, chi2, dof: %f,%f,%d"%(fit.Q, fit.chi2, fit.dof))
        print("t0 is {}\nfit-range [{},{}]".format(t0, TFIT[0], TFIT[-1]))
        sys.exit(0)

    # check fit quality by adding noise
    if NOISE:
        print('\n',30*'=','add svd, prior noise',30*'=', file=l)
        noisy_fit = fitter.lsqfit(data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT ,noise=True)
        #noisy_fit = fitter.lsqfit(data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT                                      ,add_svdnoise=False,  add_priornoise=True)

        if noisy_fit.Q < 0.1:
            print(" BAD NOISY FIT - ABORTING ")
            print("Q, chi2, dof: %f,%f,%d"%(noisy_fit.Q, noisy_fit.chi2, noisy_fit.dof))
            sys.exit(0)
        else:
            print(" GOOD FIT YOU GO BOY! ")

        print(noisy_fit.format(pstyle=None), file=l)
        dE = fit.p[tag+'dE'][:2]
        noisy_dE = noisy_fit.p[tag+'dE'][:2]
        dEo = fit.p[otag+'dE'][:2]
        noisy_dEo = noisy_fit.p[otag+'dE'][:2]
        print('      dE:', dE, file=l)
        print('noisy dE:', noisy_dE, file=l)
        print('          ', gv.fmt_chi2(gv.chi2(dE - noisy_dE)),file=l)
        print('      dEo:', dEo,file=l)
        print('noisy dEo:', noisy_dEo,file=l)
        print('          ', gv.fmt_chi2(gv.chi2(dEo - noisy_dEo)),file=l)
 
    # print results
    print_results(fit, basis, prior, data, l)

    # save masses and amplitudes and the entire lsqfit obj
    if SAVE_RES:
        save_info = {'data':data, 'prior':prior, 'svdcut':fit.correction, 'fit.p':fit.p}
        gv.dump(save_info, save_res_name, add_dependencies=True)
        with bz2.BZ2File(save_fit, 'wb') as f:
            pickle.dump( fit, f ) 


################################ FUNCTIONS ###########################################
def make_data(filename):
    dset = gv.dataset.Dataset(filename)
    print("%d cfgs"%(dset.samplesize),file=l)
    dset = gv.dataset.bin_data(dset, binsize=bin_size)
    data = c_hack*gv.dataset.avg_data(dset, warn=True)        
    return data

def make_models():
    models = []
    for i, s1 in enumerate(SRCs):
        for s2 in SRCs[i:]:
            print(i,s1,s2)
            tfit=TFIT
            if s1 == s2:
                otherdata=None
            else:
                otherdata=KEYFMT.format(s1=s2, s2=s1)
                tfit=TFIT[:D]
            if OSC:
                models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=(tag+s1,otag+s1), b=(tag+s2,otag+s2), dE=(tag+'dE',otag+'dE'), otherdata=otherdata, s=s_coeff))
            #else: models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=tag+s1, b=tag+s2, dE=tag+'dE', otherdata=otherdata))
    print("number of models = %d\n"%len(models))
    return models


def make_models_eig():
    models = []
    for i, s1 in enumerate(EIG_SRCs):
        for s2 in EIG_SRCs[i:]:
            print(i,s1,s2)
            tfit=TFIT
            if s1 == s2:
                otherdata=None
            else:
                otherdata=KEYFMT.format(s1=s2, s2=s1)
                tfit=TFIT[:D]
            if OSC:
                models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=(tag+s1,otag+s1), b=(tag+s2,otag+s2), dE=(tag+'dE',otag+'dE'), otherdata=otherdata, s=s_coeff))
            #else: models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=tag+s1, b=tag+s2, dE=tag+'dE', otherdata=otherdata))
    print("number of models = %d\n"%len(models))
    return models

def make_prior(N, basis):
    prior = basis.make_prior(nterm=N, keyfmt=(tag+'{s1}',otag+'{s1}'),  eig_srcs=EIG,)#states=([],[]))
#    prior['pp.dE'][0] = gv.gvar('3.5(0.5)')/ainv[e_str]
    return prior

def print_results(fit, basis, prior, data, logfile=None):
    print(fit.format(pstyle='m'), file=logfile)
    print(30 * '=', 'Results', 30*'=','\n', file=logfile)
    print(basis.tabulate(fit.p, keyfmt=(tag+'{s1}',otag+'{s1}'), nterm=4), file=logfile)

    E = np.cumsum(fit.p[tag+'dE'])
    Eo = np.cumsum(fit.p[otag+'dE'])
    
    outputs = collections.OrderedDict()
    outputs[tag+'grd'] = E[0]
    outputs[otag+'grd'] = Eo[0]

    inputs=collections.OrderedDict()
    inputs['prior'] = prior
    inputs['data']  = data
    inputs['svdcut'] = fit.correction

    print(gv.fmt_values(outputs), file=logfile)
    print(gv.fmt_errorbudget(outputs, inputs, colwidth=18, verify=True),file=logfile)


if __name__ == '__main__':
    gv.ranseed(1234)
    if True:
        main()
    else:
        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()
        main()
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print (s.getvalue())
