import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
from results import plym_onemm, onemm_mass, onemm_exp

hbarc = 0.197326968

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5, 4.5)))
plt.gca().tick_params(right=False, left=False, top=True,direction='in')

plt.errorbar(y=1.7, x=plym_onemm[0].mean, xerr=plym_onemm[0].sdev ,yerr=0, fmt='h',mfc='none',color='r',label=r'This work',lw=1)
plt.errorbar(y=1.6, x=onemm_mass['hadspec'].mean ,yerr=0,xerr=onemm_mass['hadspec'].sdev ,fmt='o',mfc='none',color='b',label=r'HadSpec [1610.01073]',lw=1)
plt.errorbar(y=1.5, x=onemm_mass['brambilla'].mean ,yerr=0,xerr=onemm_mass['brambilla'].sdev ,fmt='+',mfc='none',color='b',label=r'Brambilla \emph{et al.} [1805.07713]',lw=1)
plt.errorbar(y=1.4, x=onemm_mass['chen'].mean ,yerr=0,xerr=onemm_mass['chen'].sdev ,fmt='x',mfc='none',color='b',label=r'Chen \emph{et al.} [1604.03401]',lw=1)
## Experiment
plt.errorbar(y=1.2, x=onemm_exp["psi4230"].mean ,yerr=0,xerr=onemm_exp["psi4230"].sdev ,fmt='',mfc='none',color='black',lw=1)
plt.errorbar(y=1.2, x=onemm_exp["psi4360"].mean ,yerr=0,xerr=onemm_exp["psi4360"].sdev ,fmt='',mfc='none',color='black',lw=1)
plt.errorbar(y=1.2, x=onemm_exp["psi4415"].mean ,yerr=0,xerr=onemm_exp["psi4415"].sdev ,fmt='',mfc='none',color='black',lw=1)
plt.errorbar(y=1.32, x=4 ,yerr=0,xerr=1 ,fmt='-',mfc='none',color='black',lw=1.5)
##----------------------------------------------------------------------------
plt.text(4.26, 1.25, 'EXPERIMENT', fontsize=6, bbox={'facecolor': 'pink', 'alpha': 0.7, 'pad': 4})

plt.text(onemm_exp["psi4230"].mean-onemm_exp["psi4230"].sdev, 1.12, r'$\psi(4230)$', fontsize=10)
plt.text(onemm_exp["psi4360"].mean-1.3*onemm_exp["psi4360"].sdev, 1.12, r'$\psi(4360)$', fontsize=10)
plt.text(onemm_exp["psi4415"].mean-2.5*onemm_exp["psi4415"].sdev, 1.12, r'$\psi(4415)$', fontsize=10)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=10,loc='upper left')

#plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.xlabel('mass\n [GeV]', rotation=0, labelpad=0, fontsize=14)
plt.ylim(bottom=1.0, top=2.1)
plt.xlim(left=4.2, right=4.5)
frame1 = plt.gca()
frame1.axes.yaxis.set_ticklabels([])

plt.tight_layout()
plt.savefig('../../figures/onemm_fine_summary.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

