import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from commonweal import plym_onemm_gev, onemm_mass
from commonweal import a # list of lattice spacings coarsest --> finest

hbarc = 0.197326968


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5/1.618, 4.5)))
plt.gca().tick_params(right=True,top=True,direction='in')

plt.errorbar( x=plym_onemm_gev[0].mean  , y=1 , xerr=plym_onemm_gev[0].sdev ,yerr=0, fmt='h',mfc='none',color='r',label=r'This work',lw=1)

plt.errorbar(y=2, x=onemm_mass['brambilla'].mean ,yerr=0,xerr=onemm_mass['brambilla'].sdev ,fmt='o',mfc='none',color='g',label=r'Brambilla [1908.05179], EFT',lw=1)

plt.errorbar(y=3, x=onemm_mass['hadspec'].mean ,yerr=0,xerr=onemm_mass['hadspec'].sdev ,fmt='o',mfc='none',color='b',label=r'HadSpec [1610.01073], $m_\pi$=240,400 MeV',lw=1)

#---------------------------------------------------------------------------
plt.text(4.3, 1.5, 'PRELIMINARY', fontsize=6, bbox={'facecolor': 'green', 'alpha': 0.7, 'pad': 4})

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=6,loc='upper left')

#plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.xlabel('mass\n [GeV]', rotation=0, labelpad=10, fontsize=14)
plt.title(r'Summary of $\overline{c}c$ $1^{--}$ masses')
plt.ylim(bottom=0.9, top=3.5)
#plt.xlim(left=-0.0002)
frame1 = plt.gca()
frame1.axes.yaxis.set_ticklabels([])

plt.tight_layout()
plt.savefig('../../figures/onemm_fine_summary.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

