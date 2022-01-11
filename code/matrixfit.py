from __future__ import print_function   # makes this work for python2 and 3

import sys
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import os
import datetime
from parameters import corr, corrpath, l
from parameters import T, s_coeff, t0, c_hack, NEXP, TFIT, bin_size
from parameters import OSC, NOISE, WRITE_LOG
from parameters import SRCs, KEYFMT, tag, ttag, otag

# -----------------------------------------------------------------------------------
SHOWPLOTS = False         # display plots at end of fits
# ------------------------------------------------------------------------------------

print('Temporal extent: ', T)

TDATA = range(T)
tp=T

def main():
    data, basis, SVDCUT = make_data(corrpath)
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    print("svd cut is: ", SVDCUT)
    for N in NEXP:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        #print(fit.format(pstyle=None if N < 12 else 'm'))
        p0 = fit.pmean
    
    #print(prior)
    if fit.chi2/fit.dof > 1.5:
        print("reduced chi-squared: ", fit.chi2/fit.dof)
        print("t0 is {}\nfit-range [{},{}]".format(t0, TFIT[0], TFIT[-1]))
        #print_results(fit, basis, prior, data, l)
    else:
        print_results(fit, basis, prior, data, l)

    if SHOWPLOTS:
        fit.show_plots(view='ratio')

    # check fit quality by adding noise
    if NOISE:
        print('\n==================== add svd, prior noise')
        noisy_fit = fitter.lsqfit(data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT                                                   ,add_svdnoise=True	, add_priornoise=True)
        print(noisy_fit.format(pstyle=None))
        dE = fit.p[tag+'dE'][:3]
        noisy_dE = noisy_fit.p[tag+'dE'][:3]
        print('      dE:', dE)
        print('noisy dE:', noisy_dE)
        print('          ', gv.fmt_chi2(gv.chi2(dE - noisy_dE)))

def make_data(filename):
    dset = gv.dataset.Dataset(filename)
    dset = gv.dataset.bin_data(dset, binsize=bin_size)
    #data = c_hack*gv.dataset.avg_data(cf.read_dataset(filename, grep=tag))        
    data = c_hack*gv.dataset.avg_data(dset)        
    s = gv.dataset.svd_diagnosis(cf.read_dataset(filename, grep=tag), models=make_models())
    basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRCs, t=(t0, t0+1), tdata=TDATA, osc=OSC)
    return data, basis, s.svdcut

def make_models():
    models = []
    for i, s1 in enumerate(SRCs):
        for s2 in SRCs[i:]:
            tfit=TFIT if s1 == s2 else TFIT[:]
            otherdata = None if s1 == s2 else KEYFMT.format(s1=s2, s2=s1)
            if OSC:
                models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=(tag+s1,otag+s1), b=(tag+s2,otag+s2), dE=(tag+'dE',otag+'dE'), otherdata=otherdata, s=s_coeff))
            else: models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=tag+s1, b=tag+s2, dE=tag+'dE', otherdata=otherdata))
            
    return models

#def make_prior(N, basis):
#    prior = basis.make_prior(nterm=N, keyfmt=tag+'{s1}')#, states=[2]) 
    #prior['log(onemm.dE)'][2] = gv.log(gv.gvar('2.0(2)'))

    #if OSC:
        #prior1 = basis.make_prior(nterm=N, keyfmt=otag+'{s1}')#, states=[]) 
        #prior.update(prior1)    
#    return prior

def make_prior(N, basis):
   prior = basis.make_prior(nterm=N, keyfmt=tag+'{s1}', states=[])
   #prior1 = collections.OrderedDict()    
   prior['onemm.R'][0]       = gv.gvar('-0.16(1)')
   prior['onemm.r'][0]       = gv.gvar('-0.069(1)')
   #prior['onemm.R'][1]       = gv.gvar('0.19(2)')
   prior['onemm.H'][3]       = gv.gvar('-0.09(2)')
   prior['onemm.h'][3]       = gv.gvar('-0.028(3)')
   prior['log(onemm.dE)'][0] = gv.log(gv.gvar('1.41(2)'))
   #prior['log(onemm.dE)'][1] = gv.log(gv.gvar('0.50(5)'))
   #prior['log(onemm.dE)'][2] = gv.log(gv.gvar('0.30(5)'))
   
   #----------OSC PARAMS-------------#
   if OSC:
     prior['o.r'][0]       = gv.gvar('0.022(3)')
     prior['o.R'][0]       = gv.gvar('0.058(9)')
     prior['log(o.dE)'][0] = gv.log(gv.gvar('1.61(2)'))   
   return prior       
   
def print_results(fit, basis, prior, data, logfile=None):
    
    if WRITE_LOG:
      l.write('Parameters used: ' + '\n')
      l.write('file: '+ corr + '\n')
      l.write('t0: '+ str(t0) + '\n')
      l.write('T: '+ str(T) + '\n')
      l.write('tmin: '+ str(TFIT[0]) + '\n')
      l.write('tmax: '+ str(TFIT[-1]) + '\n\n')
      #l.write('svdcut: '+ str(SVDCUT) + '\n')
      #l.write('offtmin: '+ str(offtmin) + '\n')
      #l.write('offtmax: '+ str(offtmax) + '\n\n')
      l.write(30 * '=' + '\n' + 'nterm = ' +  str(NEXP[-1]) + '\n')
    print(fit.format(pstyle='m'), file=logfile)
    print(30 * '=', 'Results\n', file=logfile)
    print(basis.tabulate(fit.p, keyfmt=tag+'{s1}', nterm=6), file=logfile)
    print(basis.tabulate(fit.p, keyfmt=tag+'{s1}', nterm=6, eig_srcs=True), file=logfile)
    if OSC:
      print('Oscillating states:', file=logfile)
      print(basis.tabulate(fit.p, keyfmt=otag+'{s1}', nterm=4), file=logfile)
    E = np.cumsum(fit.p[tag+'dE'])
    
    print('Priors:', file=logfile)
    for k in [otag+SRCs[i] for i in range(len(SRCs))]:
        print('{:10}{}...'.format(k, list(prior[k][:5])), file=logfile)
    for k in [tag+SRCs[i] for i in range(len(SRCs))]:
        print('{:10}{}...'.format(k, list(prior[k][:5])), file=logfile)

    prior_eig = basis.apply(prior, keyfmt=tag+'{s1}') 
    print('\n')
    for k in [tag+str(i) for i in range(len(SRCs))]:
        print('{:10}{}...'.format(k, list(prior_eig[k][:5])), file=logfile)
    print('\n')
    if WRITE_LOG:
		    l.close()
   
        
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
