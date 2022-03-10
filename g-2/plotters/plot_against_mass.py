# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from results import vc_up_noqed, vc_down_noqed, vc_amu_up, vc_amu_down
from results import vc_up_phys_diff, vc_down_phys_diff, vc_up_diff, vc_down_diff
from results import vc_phys15
from results import vc_amu
from gvar import log
import sys

vc_amus_diff = vc_amu["ms"][2.5][1]
vc_amus      = vc_amu["ms"][2.5][0]

vc_up_phys_diff = vc_phys15['up-charge'] - vc_phys15['up-nocharge']
vc_down_phys_diff = vc_phys15['down-charge'] - vc_phys15['down-nocharge']

vc_diff = [ vc_amu["3ml"][2][1],vc_amu["5ml"][2][1],vc_amu["7ml"][2][1] ] 
vcamu = [ vc_amu["3ml"][2][0],vc_amu["5ml"][2][0],vc_amu["7ml"][2][0] ] 

#####################################################################

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':14})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
ax2.set_xticks([4/3,1,3,5,7,10])
ax2.set_xticklabels(['$m_d$','$m_l$','$3m_l$','$5m_l$','$7m_l$','$m_s$'])


ax1.set_title('0.15fm ensemble - 1844 cfgs')

ax1.errorbar(10, vc_amus_diff.mean, xerr=0, yerr=vc_amus_diff.sdev, fmt='h',color='green',lw=1 , alpha=0.8, marker='o', markersize=6)
ax2.errorbar(10, vc_amus.mean, xerr=0, yerr=vc_amus.sdev,color='green',lw=1, alpha=0.8, marker='o', markersize=6)
#ax1.errorbar(4/3, vc_down_phys_diff.mean, xerr=0, yerr=vc_down_phys_diff.sdev, fmt='h',color='green',lw=1 , alpha=0.8, marker='o', markersize=6)
ax2.errorbar(4/3, vc_down_noqed.mean, xerr=0, yerr=vc_down_noqed.sdev,color='green',lw=1, alpha=0.8, marker='o', markersize=6)

ax1.errorbar([3, 5, 7], [y.mean for y in vc_down_diff], xerr=0, yerr=[y.sdev for y in vc_down_diff], fmt='h',color='red',lw=1, alpha=0.8)
ax2.errorbar([3, 5, 7], [y.mean for y in vc_amu_down], xerr=0, yerr=[y.sdev for y in vc_amu_down], fmt='h',color='red',lw=1, alpha=0.8)
#ax2.set_xlabel(r'$\frac{m_q}{m_l}$', fontsize=25, labelpad=20)
#ax1.set_ylabel('fractional\n error\n of $\delta a_{\mu}$', rotation=0, labelpad=40)
#ax2.set_ylabel('fractional\n error\n of $a_{\mu}$', rotation=0, labelpad=40)
#ax1.yaxis.set_label_coords(-0.15,0.2)
#ax2.yaxis.set_label_coords(-0.15,0.2)
#ax1.set_yticklabels(fontsize=12)
#ax1.legend(frameon=False,fontsize=15,loc='upper right')
plt.tight_layout()
plt.savefig('../figures/vc_summary.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()
