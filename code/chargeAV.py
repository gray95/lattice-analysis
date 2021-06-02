# charge averaging on the quenched qed correlators

import gvar as gv
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

gpl_base_dir   = '../data/qed'
gpl_neu        = 'coarse/m0.01288_Pseuodoscalar_Q0.0_outcorr.gpl'
gpl_neg   		 = 'coarse/m0.01288_Pseuodoscalar_Q-0.101_outcorr.gpl'
gpl_pos        = 'coarse/m0.01288_Pseuodoscalar_Q0.101_outcorr.gpl'
gpl_av         = 'coarse/m0.01288_Pseuodoscalar_Q0.101_chargeAV_outcorr.gpl'

corrN = os.path.join(gpl_base_dir, gpl_neu)
corrM = os.path.join(gpl_base_dir, gpl_neg)
corrP = os.path.join(gpl_base_dir, gpl_pos)
corrAV = os.path.join(gpl_base_dir, gpl_av)


corrN = gv.dataset.Dataset(corrN)
corrM = gv.dataset.Dataset(corrM)
corrP = gv.dataset.Dataset(corrP)
corrAV = gv.dataset.Dataset(corrAV)
tags = corrN.keys() + corrM.keys() + corrP.keys() + corrAV.keys()
print(tags)


#for i in range(len(corrP[tags[0]])):
#  corrAV.append(0.5*(np.add(corrP[tags[0]][i], corrM[tags[1]][i])))

dneg     = gv.dataset.avg_data(corrM[tags[1]])
dpos     = gv.dataset.avg_data(corrP[tags[2]])
dneu 		 = gv.dataset.avg_data(corrN[tags[0]])
cav      = gv.dataset.avg_data(corrAV[tags[3]])

#print(cav[3].sdev)
#print(dneg[3].sdev)
#sys.exit(0)

rt_pos = np.divide(corrP[tags[2]], corrN[tags[0]])
print(rt_pos.shape)

#rt_pos = dneg/dneu
rt_cav = cav/dneu


fig, ax = plt.subplots(figsize=(13,8))
#ax.set_yscale('log')

tlim = range(22)

#ax.errorbar([0.1+i for i in tlim], [(cav[i].sdev/dneu[i].sdev) for i in tlim], xerr=0, yerr=0, fmt='bo', label='cav/neu', alpha=0.8)
#ax.errorbar(tlim, [rt_cav[i].mean for i in tlim], xerr=0, yerr=[rt_cav[i].sdev for i in tlim], fmt='ro', label='cav/neu', alpha=0.8)
ax.errorbar([0.1+i for i in tlim], [rt_pos[i].mean for i in tlim], xerr=0, yerr=np.std(rt_pos,axis=0), fmt='bx', label='pos/neu', alpha=0.8)


#plt.ylim(bottom=0.9, top=1.1)
ax.set_xlabel('time in lattice units')
#ax.set_ylabel(r'$\frac{\langle C \rangle}{\langle C \rangle}$', rotation=0, labelpad=40, fontsize=20)
ax.legend(frameon=False, fontsize=12, loc='lower left')
#ax.set_title('Ratio of pseudoscalar corrs')

#plt.savefig('./corr_ratio_test.png', dpi=500, bbox_inches='tight')
plt.show()



