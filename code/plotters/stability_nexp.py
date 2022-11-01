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

## for the fine unphysical ensemble

nexp = range(1,10)

onempf5 = [gv.gvar('2.304(19)'),gv.gvar('2.109(12)'),gv.gvar('2.147(19)')] 
onempf5.extend(6*[onempf5[-1]])

ppf5 = [gv.gvar('2.362(59)'),gv.gvar('1.897(28)'),gv.gvar('1.63(10)')] 
ppf5.extend(6*[ppf5[-1]])

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

plt.errorbar([x for x in nexp], [y.mean for y in onempf5], xerr=0, yerr=[y.sdev for y in onempf5], fmt='x',mfc='none',color='red',label='NO',lw=1)
plt.errorbar([x for x in nexp], [y.mean for y in ppf5], xerr=0, yerr=[y.sdev for y in ppf5], fmt='x',mfc='none',color='green',label='O',lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=10,loc='upper right')

plt.xlabel('N', fontsize=15)
plt.ylabel(r'$aM$', rotation=0, labelpad=13, fontsize=16)
#plt.title(r'Continuum extrapolation of $\bar{c}c$ mass')
#plt.ylim(top=5, bottom=2.5)
#plt.xlim(left=0.00)

plt.tight_layout()
plt.savefig('../../figures/stability_nexp.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

