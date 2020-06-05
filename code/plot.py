# Converts text files to numpy arrays for further processing.

import os
import sys
import numpy as np
import re
import matplotlib.pyplot as plt
import math as m 

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
          if jnow > 0 and jinc > 0 :
             tmp[i] = m.log(jnow / jinc)
          else:
             ok = False

       now = now / no_config
       inc = inc / no_config

       if (now > 0 and inc > 0) and ok   :
          meff = m.log(now/inc)
          jerr = jackknife(tmp,no_config)
       else:
          meff = -1
          jerr = 0.0

       mm[t]     = meff
       mm_err[t] = jerr

    return tt, mm, mm_err

    