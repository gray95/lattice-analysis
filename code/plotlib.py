# matplotlib plots

import numpy as np
import matplotlib.pyplot as plt
import math as m
from parameters import corrtag
import os
import sys
import re
import gvar as gv


def jackknife(jloop,nsample):
  jmean = 0.0
  for isample in  range(0,nsample) :
    jmean += jloop[isample]

  jmean  /= nsample

  jtot = 0.0
  for isample in  range(0,nsample) :
    ttt   = (jloop[isample] - jmean )
    jtot +=  ttt * ttt

  return m.sqrt( (nsample - 1.0)/(1.0*nsample) * jtot)
  
def calc_meff(corr, nt, no_config):
    tmp  = np.zeros(( no_config ))

    tmax = int(nt / 2)  - 1
    tt     = np.zeros(( tmax ))
    mm     = np.zeros(( tmax ))
    mm_err = np.zeros(( tmax ))
  
    for t in range(0,tmax):
       tt[t] = t
       now = 0.0
       inc = 0.0 
       for i in range(0, no_config  ):
          ok = True

          now = now + corr[i,t] 
          inc = inc + corr[i,t+1] 

          #  jackknife analysis
          jnow = 0
          jinc = 0
          for j in range(0, no_config  ):
            if i != j :
              jnow = jnow + corr[j,  t] 
              jinc = jinc + corr[j, t+1] 

          jnow = jnow / (no_config - 1) 
          jinc = jinc / (no_config - 1) 
          if True: #jnow > 0 and jinc > 0 :
             tmp[i] = m.log(abs(jnow) / abs(jinc))
          else:
             ok = False

       now = now / no_config
       inc = inc / no_config

       if True: #(now > 0 and inc > 0) and ok :
          meff = m.log(abs(now)/abs(inc))
          jerr = jackknife(tmp,no_config)
       else:
          meff = 0.0
          jerr = 0.0

       mm[t]     = meff
       mm_err[t] = jerr

    return tt, mm, mm_err

    
def plot_corr(c, KEY, SAVEFIG=False):   
    print ("Reading data from %s with key %s" % (c, KEY))

    corr = gv.dataset.Dataset(c)
    corr = corr.toarray()
    print(len(corr))
    corr = corr[KEY]
    print(corr.shape)
    #corr = corr * -1

    corr_mean = np.mean(corr, axis=0)
    corr_uncertainty = np.std(corr, axis=0)
    
    no_config = corr.shape[0]
    nt        = corr.shape[1]
    
    tt, meff, meff_err = calc_meff(corr, nt, no_config) 
    
    ### plot the data
    
    plt.subplot(211)
    plt.errorbar(range(len(corr_mean)),[abs(corr_mean[q]) for q in range(len(corr_mean))], fmt='bo', yerr=None  , alpha=0.5, markersize=2)
    plt.yscale('log', nonposy='clip')
    plt.title(r'hyb-hyb $1^{--}$ mean correlator')  
    plt.ylabel(r'$\langle 0|e^{-HT}|0 \rangle$')
    
    E = []
    for i in range(len(corr_mean)-1):
        hldr = m.log(abs((corr_mean[i]/corr_mean[i+1])))
        E.append(hldr)
        
    #plt.subplot(312)
    #plt.plot(range(len(E)), [E[x] for x in range(len(E))], 'b+')
    #plt.ylabel(r'$\log( \frac{G(t)}{G(t+a)})$')
    #plt.xlabel('time')
    #plt.axis([0,nt,-6,6])
    
    plt.subplot(212)
    plt.errorbar(tt, meff , meff_err  ,  fmt= 'ro', markersize=2, alpha=0.5)
    plt.xlabel('t')
    plt.ylabel('meff')
    plt.xlim(0,nt//4)
    plt.ylim(0,6)
    plt.yticks(np.arange(0, 4, 0.5))
    plt.xticks(np.arange(0, nt//4, 1))
    plt.grid(axis='y')
    
    if SAVEFIG :
      plt.savefig(('../figures/'+corrtag+'_'+KEY+'.png'), dpi=600)
    
    
    plt.show()
        
        
