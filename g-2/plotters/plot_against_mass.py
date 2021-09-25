# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from commonweal import bmw_diff, giu_diff, bmw_diffsib
from commonweal import c_amu, c_rt, c_diff, vc_amu, vc_rt, vc_diff
from commonweal import vc_amu_phys, vc_rt_phys, vc_diff_phys

fig, ax = plt.subplots(figsize=(10,7))
ax.set_yscale('linear')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

ax.errorbar(0.95, bmw_diffsib.mean, xerr=0, yerr=bmw_diffsib.sdev, fmt='h',mfc='none',color='pink',lw=2, label=r'BMW ($m_u \neq m_d$)', alpha=0.8)
ax.errorbar(0.95, bmw_diff.mean, xerr=0, yerr=bmw_diff.sdev, fmt='h',mfc='none',color='purple',lw=2, label='BMW', alpha=0.8)
ax.errorbar(0.95, giu_diff.mean, xerr=0, yerr=giu_diff.sdev, fmt='h',mfc='none',color='orange',lw=2, label='Giusti et al.', alpha=0.8)
ax.errorbar(1.05, vc_diff_phys.mean, xerr=0, yerr=vc_diff_phys.sdev, fmt='h',mfc='none',color='green',lw=1, label=r'very coarse physical point ($m_u \neq m_d$)', alpha=0.8)

ax.errorbar([3, 5, 7], [x.mean for x in vc_diff], xerr=0, yerr=[x.sdev for x in vc_diff], fmt='h',mfc='none',color='red',lw=1, label='very coarse-185 configs', alpha=0.8)

ax.errorbar([3.05, 5.05, 7.05], [x.mean for x in c_diff], xerr=0, yerr=[x.sdev for x in c_diff], fmt='h',mfc='none',color='blue',lw=1.5, label='coarse-56 configs', alpha=0.8)

ax.set_xlabel(r'$\frac{m_q}{m_l}$', fontsize=25, labelpad=20)
ax.set_ylabel(r'$\frac{\delta a_{\mu}^{(\mathrm{light})}}{a_{\mu}^{(\mathrm{light})}}$', fontsize=26, rotation=0, labelpad=40)
plt.yticks(fontsize=12)#ax.set_yticklabels(fontsize=12)
ax.legend(frameon=False,fontsize=15,loc='upper right')
#plt.ylim(top=6)
#plt.xlim(left=1)
#plt.savefig('../figures/g-2_diff.png', dpi=500, bbox_inches="tight")
plt.tight_layout()
plt.show()

plt.close()
