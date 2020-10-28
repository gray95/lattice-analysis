from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
import corrfitter as cf
import os
import datetime

# -----------------------------------------------------------------------------------
SHOWPLOTS = False         # display plots at end of fits
OSC = True                # include oscillating states 
NOISE = False             # check fit quality by adding noise
WRITE_LOG = False					# write out a log file with results and parameters.
# ------------------------------------------------------------------------------------



#SRC = ['l', 'g']            # labels for the sources     
SRC = ['H', 'R', 'h', 'r']
KEYFMT = 'onemm.{s1}{s2}'   # keys
T = 96                      # temporal extent of lattice
TDATA = range(T)
SVDCUT = 0.0005
NEXP = range(1,13)            # number of exponentials in fit

tmin = 4               # start fit from here for diagonal elements (ll, gg,..)
tmax = 14
offtmin = 12					# off-diagonal elements (lg, gl,..)
offtmax = 18            
t0 = 4                      # initial timeslice to generate priors

c_hack =  1 									# sometimes -1 needed to generate priors

tag = KEYFMT[:-8]           # with a . "onemp."
ttag = tag[:-1]             # no .     "onemp"
otag = ttag + '_o.'         # oscillating tags "onemp_o."


corr = 'comb4_onemmHy_vec_m0.450.txt'  
 
log_folder = 'l3296f211b630m0074m037m440-coul-v5'
file = os.path.join('../data/proc_corrs', log_folder, corr)   # path to file


log_name =  'LOG' + '_' + corr       # name of log file.
if WRITE_LOG:
    l = open('../log/'+log_folder+'/'+log_name, 'w+')


def main():
    data, basis = make_data(file)
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in NEXP:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        print(fit.format(pstyle=None if N < 12 else 'm'))
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
    data = c_hack*gv.dataset.avg_data(cf.read_dataset(filename, grep=ttag))        
    basis = cf.EigenBasis(data, keyfmt=KEYFMT, srcs=SRC, t=(t0, t0+2), tdata=TDATA)
    return data, basis

def make_models():
    models = []
    for i, s1 in enumerate(SRC):
        for s2 in SRC[i:]:
            tfit=TDATA[tmin:tmax] if s1 == s2 else TDATA[offtmin:offtmax]
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
    l.write('Parameters used: ' + '\n')
    l.write('file: '+ file + '\n')
    l.write('t0: '+ str(t0) + '\n')
    l.write('T: '+ str(T) + '\n')
    l.write('tmin: '+ str(tmin) + '\n')
    l.write('tmax: '+ str(tmax) + '\n\n')
    l.write('offtmin: '+ str(offtmin) + '\n')
    l.write('offtmax: '+ str(offtmax) + '\n\n')

    l.write(30 * '=' + '\n' + 'nterm = ' +  str(N[-1]) + '\n')
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
