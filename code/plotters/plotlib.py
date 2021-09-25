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

    tmax   = nt//4  - 1
    tt     = np.zeros(( tmax ))
    mm     = np.zeros(( tmax ))
    sm_mm  = np.zeros(( tmax ))
    mm_err = np.zeros(( tmax ))
    sm_mm_err=np.zeros(( tmax ))
  
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

    for t in range(1, tmax-1):
       sm_mm[t] = 0.25*(mm[t+1]+2*mm[t]+mm[t-1])
    sm_mm[0] = mm[0]

    return tt, mm, mm_err, sm_mm

    
def plot_corr(c, SAVEFIG=False):   
    print ("Reading data from %s" % c)

    corr = gv.dataset.Dataset(c)
    keys = corr.keys()
    keys = list(keys)
    corr = corr.toarray()
    print(keys)
    KEY=keys[0]
#    corr_vt = corr[keys[4]]
    corr_hy = corr[KEY]
#    print(corr_vt.shape)
    print(corr_hy.shape)
   
    no_config = corr_hy.shape[0]
    nt        = corr_hy.shape[1]
    
#    tt_vt, meff_vt, meff_err_vt, smooth_meff_vt = calc_meff(corr_vt, nt, no_config) 
    tt_hy, meff_hy, meff_err_hy, smooth_meff_hy = calc_meff(corr_hy, nt, no_config) 
    
    ### plot the data
    
 
    #plt.subplot(212)
#    plt.errorbar(tt_vt, meff_vt, yerr=meff_err_vt, fmt= 'ro', markersize=5, alpha=0.5)
    plt.errorbar(tt_hy, meff_hy, yerr=meff_err_hy, fmt= 'ro', markersize=5, alpha=0.7)
    plt.fill_between(tt_hy, 1.993, 2.079, color='green', alpha=0.4) 
    plt.xlabel('t/a', fontsize=14)
    plt.ylabel('m_eff', rotation=0, labelpad=15, fontsize=14)
    plt.xlim(0,nt//8)
    plt.ylim(1,4)
    plt.xticks(np.arange(0, nt//8, 1))
    plt.grid(axis='x')
    plt.title('effective mass of the unsmeared $1^{-+}$ hybrid correlator')
    
    if SAVEFIG :
      plt.savefig(('../figures/'+corrtag+'_'+KEY+'.png'), dpi=600)
    
    
    plt.show()
        
        
