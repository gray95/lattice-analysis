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

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

for E in ratio:
    x = (hbarc/ainv[E])**2
    y = ratio[E]
    plt.errorbar(x.mean, y.mean ,xerr=x.sdev, yerr=y.sdev,fmt='o',mfc='none',color='r',lw=1)

#plt.errorbar([i.mean**2 for i in aa], [x.mean for x in pp_m],xerr=0,yerr=[y.sdev for y in pp_m],fmt='x',mfc='none',color='g',label='O - $1^{++}$ parity partner',lw=1)

#handles,labels = plt.gca().get_legend_handles_labels()
#handles = [h[0] for h in handles]
#plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.ylabel(r'$\frac{M_{\rm{NO}}}{M_{\rm{O}}}$', rotation=0, labelpad=15, fontsize=14)
#plt.title(r'Continuum extrapolation of $\bar{c}c$ mass')
#plt.ylim(top=5, bottom=2.5)
#plt.xlim(left=0.00)


plt.tight_layout()
plt.savefig('../../figures/onemp_ratio.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

