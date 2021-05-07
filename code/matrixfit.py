from __future__ import print_function   # makes this work for python2 and 3

import sys
import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import os
import datetime
from parameters import corr, corrpath, l
from parameters import T, s_coeff
from parameters import OSC, NOISE, WRITE_LOG
from parameters import KEYFMT, tag, ttag, otag

# -----------------------------------------------------------------------------------
SHOWPLOTS = False         # display plots at end of fits
# ------------------------------------------------------------------------------------

print('Temporal extent: ', T)

SRCs = ['l', 'g']            # labels for the sources     
TDATA = range(T)
SVDCUT = 0.0005
NEXP = range(1,13)            # number of exponentials in fit
tp=T
tmin = 3               # start fit from here for diagonal elements (ll, gg,..)
tmax = tmin+8 
offtmin = tmin					# off-diagonal elements (lg, gl,..)
offtmax = offtmin+8            
t0 = 2                      # initial timeslice to generate priors

c_hack = -1									# sometimes -1 needed to generate priors


def main():
    data, basis = make_data(corrpath)
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in NEXP:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        print(fit.format(pstyle=None if N < 12 else 'm'))
        p0 = fit.pmean
    print_results(fit, basis, prior, data, l)
    if SHOWPLOTS:
        fit.show_plots(view='ratio')

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
    dset = gv.dataset.Dataset(filename)
    data = c_hack*gv.dataset.avg_data(cf.read_dataset(filename, grep=ttag))        
    if OSC:
      basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRCs, t=(t0, t0+2), tdata=TDATA)
    else:
      basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRCs, t=(t0, t0+1), tdata=TDATA)
    #T = data[dset.keys()[0]].size
    return data, basis

def make_models():
    models = []
    for i, s1 in enumerate(SRCs):
        for s2 in SRCs[i:]:
            tfit=TDATA[tmin:tmax] if s1 == s2 else TDATA[offtmin:offtmax]
            otherdata = None if s1 == s2 else KEYFMT.format(s1=s2, s2=s1)
            if OSC:
                models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=(tag+s1,otag+s1), b=(tag+s2,otag+s2), dE=(tag+'dE',otag+'dE'), otherdata=otherdata, s=s_coeff))
            else: models.append(cf.Corr2(datatag=KEYFMT.format(s1=s1, s2=s2),tdata=TDATA, tfit=tfit, tp=tp, a=tag+s1, b=tag+s2, dE=tag+'dE', otherdata=otherdata))
            
    return models

def make_prior(N, basis):
    prior = basis.make_prior(nterm=N, keyfmt=tag+'{s1}')#, eig_srcs=True) 
    if OSC: 
        prior1 = collections.OrderedDict()
        prior1['onemp_o.l'] = gv.gvar(['1(0.9)'] + (N-1) * ['0.5(5)'])  
        prior1['onemp_o.g'] = gv.gvar(['1(0.9)'] + (N-1) * ['0.5(5)'])   
        prior1['log(onemp_o.dE)'] = gv.log(gv.gvar(['1.6(2)'] + (N-1)*['0.5(4)']))
        prior.update(prior1)    
    return prior
            
    
def print_results(fit, basis, prior, data, logfile=None):
    
    if WRITE_LOG:
      l.write('Parameters used: ' + '\n')
      l.write('file: '+ corr + '\n')
      l.write('t0: '+ str(t0) + '\n')
      l.write('T: '+ str(T) + '\n')
      l.write('tmin: '+ str(tmin) + '\n')
      l.write('tmax: '+ str(tmax) + '\n\n')
      l.write('offtmin: '+ str(offtmin) + '\n')
      l.write('offtmax: '+ str(offtmax) + '\n\n')
      l.write(30 * '=' + '\n' + 'nterm = ' +  str(NEXP[-1]) + '\n')
    print(fit.format(pstyle='m'), file=logfile)
    print(30 * '=', 'Results\n', file=logfile)
    print(basis.tabulate(fit.p, keyfmt=tag+'{s1}', nterm=6), file=logfile)
    print(basis.tabulate(fit.p, keyfmt=tag+'{s1}', nterm=6, eig_srcs=True), file=logfile)
    if OSC:
      print('Oscillating states:', file=logfile)
      print(basis.tabulate(fit.p, keyfmt=otag+'{s1}', nterm=3), file=logfile)
    E = np.cumsum(fit.p[tag+'dE'])
    #outputs = collections.OrderedDict()
    #outputs['a*E(2s-1s)'] = E[1] - E[0]
    #outputs['a*E(3s-1s)'] = E[2] - E[0]
    #outputs['E(3s-1s)/E(2s-1s)'] = (E[2] - E[0]) / (E[1] - E[0])
    #inputs = collections.OrderedDict()
    #inputs['prior'] = prior
    #inputs['data'] = data
    #inputs['svdcut'] = fit.svdcorrection
    #print(gv.fmt_values(outputs), file=logfile)
    #print(gv.fmt_errorbudget(outputs, inputs, colwidth=18), file=logfile)

    print('Prior:\n', file=logfile)
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
