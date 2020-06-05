# matplotlib plots

import numpy as np
from plot import jackknife, calc_meff
import matplotlib.pyplot as plt
import math as m

from parameters import corrtag

def plot_corr(c, SAVEFIG):   
    print ("Reading data from  " , c)
    data = open(c, 'r')
    
    corr = np.genfromtxt(data)
    #corr = corr * -1
    
    corr_mean = np.mean(corr, axis=0)
    uncertainty = np.std(corr, axis=0)
    
    corr_mean = corr_mean[1:]
    uncertainty = uncertainty[1:]
    
    no_config = corr.shape[0]
    nt        = corr.shape[1]
    
    tt, meff, meff_err = calc_meff(corr, nt, no_config) 
    
    ### plot the data
    
    plt.subplot(211)
    plt.errorbar(range(len(corr_mean)),[abs(corr_mean[q]) for q in range(len(corr_mean))], fmt='bo', yerr=None  , alpha=0.5, markersize=2)
    plt.yscale('log', nonposy='clip')
    plt.title(r'$1^{-+}$ mean correlator')  
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
    plt.ylim(-2,6)
    plt.yticks(np.arange(-2, 6, 1))
    plt.xticks(np.arange(0, nt//4, 1))
    plt.grid(axis='y')
    
    if SAVEFIG :
      plt.savefig((corrtag + '.png'), dpi=600)
    
    
    plt.show()
        
        