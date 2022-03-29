# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from common_results import bmw_amu_s, hpqcd_amu_s, vc_amu_s, c_amu_s
from common_results import c_amudiffrt_s, vc_amudiffrt_s, giu_amudiffrt_s

fig, ax = plt.subplots(figsize=(5,7))
ax.set_yscale('linear')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

ax.errorbar(vc_amudiffrt_s.mean, 0.3, xerr=vc_amudiffrt_s.sdev, yerr=0, fmt='h',mfc='none',color='red',lw=3, label='185 configs 0.15fm', alpha=0.8)
ax.errorbar(c_amudiffrt_s.mean, 0.6, xerr=c_amudiffrt_s.sdev, yerr=0, fmt='h',mfc='none',color='orange',lw=3, label='56 configs 0.12fm', alpha=0.8)
#ax.errorbar(hpqcd_amu_s.mean, 1, xerr=hpqcd_amu_s.sdev, yerr=0, fmt='h',mfc='none',color='green',lw=3, label='HPQCD', alpha=0.8)
#ax.errorbar(bmw_amu_s.mean, 1.4, xerr=bmw_amu_s.sdev, yerr=0, fmt='h',mfc='none',color='black',lw=3, label='BMW', alpha=0.8)
#ax.errorbar(giu_amudiffrt_s.mean, 0.6, xerr=giu_amudiffrt_s.sdev, yerr=0, fmt='h',mfc='none',color='orange',lw=3, label='Giusti et al.', alpha=0.8)


ax.set_yticklabels([])
ax.set_xlabel(r'$\frac{\delta a_{\mu}^{(s)}}{a_{\mu}^{(s)}}$', fontsize=12, labelpad=-5)
#ax.set_ylabel('$a_{\mu}^{[qcd+qed]}$\n\nx$10^{-8}$', fontsize=22, rotation=0, labelpad=50)
ax.legend(frameon=False,fontsize=15,loc='upper right')
plt.ylim(top=0.8)
#plt.xlim(right=0)
#plt.savefig('../figures/g-2_diffstrange.png', dpi=500, bbox_inches="tight")
#plt.tight_layout()
plt.show()

plt.close()
