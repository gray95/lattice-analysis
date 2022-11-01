import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from results import ratio
from commonweal import ainv, hbarc
from parameters import save_res_name
import collections

## plot value of quantity vs. number of exponentials

## for the coarse physical ensemble

nexp = range(1,10)

chi2_onempc = [320,3.1,0.66,0.58,0.53,0.52,0.52,0.52,0.53] 
Q_onempc    = [0,9.1e-05,0.82,0.89,0.92,0.93,0.92,0.92,0.92]

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

plt.plot([x for x in nexp[1:]], [y for y in chi2_onempc[1:]],'bx',label=r'$\chi^2_\nu$',lw=1)
plt.plot([x for x in nexp[1:]], [y for y in Q_onempc[1:]],'ro',label=r'$Q$',lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
#handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=10,loc='upper right')

plt.xlabel('N', fontsize=15)
#plt.ylabel(r'$\chi^2$', rotation=0, labelpad=10, fontsize=16)
#plt.title(r'Continuum extrapolation of $\bar{c}c$ mass')
#plt.ylim(top=5, bottom=2.5)
#plt.xlim(left=0.00)

plt.tight_layout()
plt.savefig('../../figures/stability_nexp_chi2.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

