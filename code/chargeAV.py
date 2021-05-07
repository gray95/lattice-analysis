# charge averaging on the quenched qed correlators

import gvar as gv
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

gpl_base_dir = '../data/qed'
gpl_pos         = 'vcoarse/dan_b/m0.01698_Q0.202_G5G5_notsm.gpl'
gpl_neg         = 'vcoarse/dan_b/m0.01698_Q-0.202_G5G5_notsm.gpl'
gpl_neu   			= 'vcoarse/dan_b/m0.01698_Q0.0_G5G5_notsm.gpl'

#gpl_pos         = 'vcoarse/dan_b/m0.007278_Q0.101_GXGX_notsm.gpl'
#gpl_neg         = 'vcoarse/dan_b/m0.007278_Q-0.101_GXGX_notsm.gpl'
#gpl_neu   			= 'vcoarse/dan_b/m0.007278_Q0.0_GXGX_notsm.gpl'

corrP = os.path.join(gpl_base_dir, gpl_pos)
corrM = os.path.join(gpl_base_dir, gpl_neg)
corrN = os.path.join(gpl_base_dir, gpl_neu)

corrP = gv.dataset.Dataset(corrP)
corrM = gv.dataset.Dataset(corrM)
corrN = gv.dataset.Dataset(corrN)
corrAV = []
tags = corrP.keys() + corrM.keys()  + corrN.keys() #+ corr3.keys()

for i in range(len(corrP[tags[0]])):
  corrAV.append(0.5*(np.add(corrP[tags[0]][i], corrM[tags[1]][i])))

#cav = 0.5*np.add(corrP[tags[0]], corrM[tags[1]])
rt_pos = np.divide(corrP[tags[0]], corrN[tags[2]])
rt_cav = np.divide(corrAV, corrN[tags[2]])

dneg     = gv.dataset.avg_data(corrM[tags[1]])
dpos     = gv.dataset.avg_data(corrP[tags[0]])
dneu 		 = gv.dataset.avg_data(corrN[tags[2]])
cav      = gv.dataset.avg_data(corrAV)

print(dpos[:10])
print(dneg[:10])
print(cav[:10])
print()


rt_cav = [(cav[i]/dneu[i]) for i in range(len((cav)))]
rt_pos = [(dpos[i]/dneu[i]) for i in range(len((dpos)))]
#rt_pos = dpos/dneu
#rt_cav = cav/dneu

#sys.exit(0)

fig, ax = plt.subplots(figsize=(13,8))
#ax.set_yscale('log')

tlim = range(20)

ax.errorbar([0.1+i for i in tlim], [rt_pos[i].mean for i in tlim], xerr=0, yerr=rt_pos[i].sdev, fmt='bo', label='P', alpha=0.8)
ax.errorbar(tlim, [rt_cav[i].mean for i in tlim], xerr=0, yerr=rt_cav[i].sdev, fmt='ro', label='AV', alpha=0.8)
#ax.errorbar([+0.1+i for i in tlim], [dneu[i].mean for i in tlim], xerr=0, yerr=0, fmt='rx', label='neu', alpha=0.8)


#plt.ylim(bottom=0.9, top=1.1)
ax.set_xlabel('time in lattice units')
#ax.set_ylabel(r'$\frac{\langle C \rangle}{\langle C \rangle}$', rotation=0, labelpad=40, fontsize=20)
ax.legend(frameon=False, fontsize=12, loc='lower left')
#ax.set_title('Ratio of pseudoscalar corrs')

#plt.savefig('../g-2/figures/cav_corr_ratio_test.png', dpi=500, bbox_inches='tight')
plt.show()



