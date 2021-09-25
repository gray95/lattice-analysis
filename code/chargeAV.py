# charge averaging on the quenched qed correlators

import gvar as gv
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

gpl_base_dir   = '../data/qqed'
gpl_neu        = 'vcoarse/b_full/m0.0677_Pseuodoscalar_Q0.0_outcorr.gpl'
gpl_neg   		 = 'vcoarse/b_full/m0.0677_Pseuodoscalar_Q-0.202_outcorr.gpl'
gpl_pos        = 'vcoarse/b_full/m0.0677_Pseuodoscalar_Q0.202_outcorr.gpl'
gpl_av         = 'vcoarse/b_full/m0.0677_Pseuodoscalar_Q0.202_chargeAV_outcorr.gpl'

corrN = os.path.join(gpl_base_dir, gpl_neu)
corrM = os.path.join(gpl_base_dir, gpl_neg)
corrP = os.path.join(gpl_base_dir, gpl_pos)
corrAV = os.path.join(gpl_base_dir, gpl_av)

corrN = gv.dataset.Dataset(corrN)
corrM = gv.dataset.Dataset(corrM)
corrP = gv.dataset.Dataset(corrP)
corrAV = gv.dataset.Dataset(corrAV)

NKeys = corrN.keys() 
MKeys = corrM.keys() 
PKeys = corrP.keys() 
AVKeys= corrAV.keys()

CAV=[]
for j in range(len(corrM)):
	tmp=np.average( (np.array( [corrP[PKeys[0]], corrM[MKeys[0]]] )), axis=0)
	print(tmp.shape)
	CAV.append(tmp)	

dneg     = gv.dataset.avg_data(corrM[MKeys[0]])
dpos     = gv.dataset.avg_data(corrP[PKeys[0]])
dneu 		 = gv.dataset.avg_data(corrN[NKeys[0]])
mycav    = gv.dataset.avg_data(CAV[0])
cav      = gv.dataset.avg_data(corrAV[AVKeys[0]])

#print(cav[3].sdev)
#print(dneg[3].sdev)

rt_cav = cav/mycav
rt_pos = dpos/dneu

fig, ax = plt.subplots(figsize=(13,8))
#ax.set_yscale('log')

tlim = range(12)


ax.errorbar([-0.1+i for i in tlim], [rt_cav[i].mean for i in tlim], xerr=0, yerr=[rt_cav[i].sdev for i in tlim], fmt='ro', label='averaged/neu', alpha=0.8)
ax.errorbar([0.1+i for i in tlim], [rt_pos[i].mean for i in tlim], xerr=0, yerr=[rt_pos[i].sdev for i in tlim], fmt='bx', label='pos/neu', alpha=0.8)

#plt.ylim(bottom=0.9, top=1.1)
ax.set_xlabel('time in lattice units')
#ax.set_ylabel(r'$\frac{\langle C \rangle}{\langle C \rangle}$', rotation=0, labelpad=40, fontsize=20)
ax.legend(frameon=False, fontsize=12, loc='lower left')
#ax.set_title('Ratio of pseudoscalar corrs')

#plt.savefig('./corr_ratio_test.png', dpi=500, bbox_inches='tight')
plt.show()



