from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import matplotlib
import os
import datetime

# -----------------------------------------------------------------------------------
SHOWPLOTS = False         # display plots at end of fits
OSC = True                # include oscillating states 
NOISE = False             # check fit quality by adding noise
WRITE_LOG = True          # write out a log file with results and parameters.
# ------------------------------------------------------------------------------------


SRC = ['l', 'g']            # labels for the sources     
KEYFMT = 'onemp.{s1}{s2}'   # keys
T = 96                      # temporal extent of lattice
TDATA = range(T)
SVDCUT = 0.001
NEXP = 13                   # number of exponentials in fit
diag = 12                   # start fit from here for diagonal elements (ll, gg)
offdiag = 28                # fit for off diagonal elements (gl, lg etc.)
t0 = 9                      # initial timeslice to generate priors

tag = KEYFMT[:-8]           # with a . "onemp."
ttag = tag[:-1]             # no .     "onemp"
otag = ttag + '_o.'         # oscillating tags "onemp_o."

corr = 'Hy_m0.450_doublysmeared_lrho.txt'                                       # filename
file = os.path.join('proc_corrs', 'l3296f211b630m0074m037m440-coul-v5', corr)   # path to file


log_name =  'LOG' + '_' + corr       # name of log file.
if WRITE_LOG:
    l = open('log/'+log_name, 'w+')


def main():
    data, basis = make_data(file)
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in [NEXP]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        print(fit.format(pstyle=None if N < 7 else 'm'))
        p0 = fit.pmean
    print_results(fit, basis, prior, data)
    if WRITE_LOG:
        write_results(fit, basis, prior, data, NEXP)
    if SHOWPLOTS:
        fit.show_plots(save='etac.{}.png', view='ratio')

    # check fit quality by adding noise
    #print('\n==================== add svd, prior noise')
    #noisy_fit = fitter.lsqfit(
    #    data=data, prior=prior, p0=fit.pmean, svdcut=None,
    #    add_svdnoise=True, add_priornoise=True,
    #   )
    #print(noisy_fit.format(pstyle=None))
    #dE = fit.p['etab.dE'][:3]
    #noisy_dE = noisy_fit.p['etab.dE'][:3]
    #print('      dE:', dE)
    #print('noisy dE:', noisy_dE)
    #print('          ', gv.fmt_chi2(gv.chi2(dE - noisy_dE)))

def make_data(filename):
    data = gv.dataset.avg_data(cf.read_dataset(filename, grep=tag[:-1]))            # sometimes *-1 forces fit
    basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRC, t=(t0, t0+2), tdata=TDATA)
    return data, basis

def make_models():
    models = []
    for i, s1 in enumerate(SRC):
        for s2 in SRC[i:]:
            tfit = TDATA[diag:-diag] if s1 == s2 else TDATA[offdiag:-offdiag]
            otherdata = None if s1 == s2 else KEYFMT.format(s1=s2, s2=s1)
            if OSC:
                models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=T,
                    a=(tag+s1,otag+s1), b=(tag+s2,otag+s2), dE=(tag+'dE',otag+'dE'),
                    otherdata=otherdata))
            else: models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=T,
                    a=tag+s1, b=tag+s2, dE=tag+'dE', otherdata=otherdata))
            
    return models

def make_prior(N, basis):
    prior = basis.make_prior(nterm=N, keyfmt=tag+'{s1}') 
    if OSC: 
        prior1 = basis.make_prior(nterm=N, keyfmt=otag+'{s1}')  
        prior.update(prior1)    
    return prior
            
    
def print_results(fit, basis, prior, data):
    print(30 * '=', 'Results\n')
    print(basis.tabulate(fit.p, keyfmt=tag+'{s1}'))
    if OSC:
      print(basis.tabulate(fit.p, keyfmt=otag+'{s1}'))
    #print(basis.tabulate(fit.p, keyfmt=tag+'{s1}', eig_srcs=True))
    E = np.cumsum(fit.p[tag+'dE'])
    outputs = collections.OrderedDict()
    outputs['a*E(2s-1s)'] = E[1] - E[0]
    outputs['a*E(3s-1s)'] = E[2] - E[0]
    outputs['E(3s-1s)/E(2s-1s)'] = (E[2] - E[0]) / (E[1] - E[0])
    inputs = collections.OrderedDict()
    inputs['prior'] = prior
    inputs['data'] = data
    inputs['svdcut'] = fit.svdcorrection
    print(gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs, colwidth=18))
    print('Prior:\n')
    for k in [tag+SRC[0], tag+SRC[1]]:
        print('{:13}{}'.format(k, list(prior[k])))
    print()
    #prior_eig = basis.apply(prior, keyfmt='onemp.{s1}')
    #for k in [tag+SRC[0], tag+SRC[1]]:
    #    print('{:13}{}'.format(k, list(prior_eig[k])))
    
def write_results(fit, basis, prior, data, N):
    l.write(30 * '=' + '\n' + 'nterm = ' +  str(N) + '\n')
    l.write(fit.format(pstyle=None if N < 7 else 'm'))
    
    l.write(30 * '=' + ' Results\n')
    l.write(basis.tabulate(fit.p, keyfmt=tag+'{s1}'))
    if OSC:
      l.write(basis.tabulate(fit.p, keyfmt=otag+'{s1}'))
    E = np.cumsum(fit.p[tag+'dE'])
    outputs = collections.OrderedDict()
    outputs['a*E(2s-1s)'] = E[1] - E[0]
    outputs['a*E(3s-1s)'] = E[2] - E[0]
    outputs['E(3s-1s)/E(2s-1s)'] = (E[2] - E[0]) / (E[1] - E[0])
    inputs = collections.OrderedDict()
    inputs['prior'] = prior
    inputs['data'] = data
    inputs['svdcut'] = fit.svdcorrection
    l.write(gv.fmt_values(outputs))
    l.write(gv.fmt_errorbudget(outputs, inputs, colwidth=18))
    l.write('Prior:\n')
    for k in [tag+SRC[0], tag+SRC[1]]:
        l.write('{:8}{}\n'.format(k, list(prior[k])))
    
        
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