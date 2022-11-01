import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from results import onemp_m, pp_m
from commonweal import ainv, hbarc

a = [ 1/ainv[x] for x in onemp_m ]
aa = [ 1/ainv[x] for x in pp_m ]
onemp_m = [ onemp_m[y] for y in onemp_m ]
pp_m = [ pp_m[y] for y in pp_m ]


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

plt.errorbar([(i.mean*hbarc)**2 for i in a], [x.mean for x in onemp_m],xerr=0,yerr=[y.sdev for y in onemp_m],fmt='h',mfc='none',color='r',label=r'This work',lw=1)
#plt.errorbar([i.mean**2 for i in a], [x.mean for x in pp_m],xerr=0,yerr=[y.sdev for y in pp_m],fmt='h',mfc='none',color='g',label='$1^{++}$ $\overline{c}c$ (parity partner)',lw=1)
plt.errorbar(0.0, y=4.53 ,xerr=0,yerr=0.15 ,fmt='h',color='r',lw=1)

# other groups
##plt.errorbar([0.12**2], [1.233+2.984],xerr=0,yerr=[16/1000],fmt='s',mfc='none',color='b',label=r'HadSpec 1204.5425, $m_\pi$=400 MeV',lw=1)
plt.errorbar([0.12**2], [1.233+2.984],xerr=0,yerr=[16/1000],fmt='o',mfc='none',color='b',lw=1)

plt.errorbar([0.12**2], [1.326+2.984],xerr=0,yerr=[23/1000],fmt='o',mfc='none',color='b',label=r'HadSpec 1204.5425,1610.01073',lw=1)

plt.errorbar([0.1145**2], [4.154],xerr=0,yerr=[54/1000],fmt='s',mfc='none',color='g',label=r'Bali et al. 1110.2381',lw=1)

plt.errorbar([0.18**2], [4.309],xerr=0,yerr=[2/1000],fmt='D',mfc='none',color='k',label=r'Ma et al. 1910.09819',lw=1)

###  decay thresholds
plt.plot([0.0 , 0.010 ] , [ 2*1.865 , 2*1.865  ] , "--", color = "black" , lw=0.4)
plt.text(0.001 , 2*1.865 +0.05  , r"D$\overline{D}$", fontsize=6, alpha=0.5)

plt.plot([0.0 , 0.025 ] , [ (2.343+2.007), (2.343+2.007) ] , "--", color = "black", lw=0.4)
plt.text(0.001 , (2.343+2.007)+0.05  , r"$D_0^{\star}(2300)\overline{D^{\star}}$", fontsize=6, alpha=0.5)

plt.plot([0.0 , 0.025 ] , [ (2.422+1.865), (2.422+1.865) ] , "--", color = "black",lw=0.4)
plt.text(0.001 , (2.422+1.865)-0.15  , r"$D_1\overline{D}$", fontsize=6, alpha=0.5)

# two P-wave 
plt.text(0.001 , 2*2.343 +0.05  , r"$D_0^{\star}(2300)\overline{D_0}^{\star}(2300)$", fontsize=6, alpha=0.5)
plt.plot([0.0 , 0.025 ] , [ 2*2.343 , 2*2.343  ] , "--", color = "black", lw=0.4 )
#---------------------------------------------------------------------------

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=10)
plt.ylabel('mass\n [GeV]', rotation=0, labelpad=10, fontsize=10)
#plt.title(r'Summary of $\overline{c}c$ $1^{-+}$ masses')
plt.ylim(top=5, bottom=3.3)
plt.xlim(left=-0.0002)


plt.tight_layout()
plt.savefig('../../figures/onemp_summary.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

