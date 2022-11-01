# pol averaging of any correlators

import gvar as gv
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


#l3248f211b580m002426m06730m8447
#l3296f211b630m0074m037m440-coul-v5
#l3264f211b600m00507m0507m628a-coul-v5
#l4864f211b600m001907m05252m6382 - coarse phys
# coarse 0.0527 0.01288 0.0092 0.00552
#mass = 'm0.0362'
basedir = '../data/hybrid'
ensemble = 'l6496f211b630m0012m0363m432'

#corrx = mass+'_Rho_GXGX_Q0.202_chargeAV_outcorr.gpl'
#corry = mass+'_Rho_GYGY_Q0.202_chargeAV_outcorr.gpl'
#corrz = mass+'_Rho_GZGZ_Q0.202_chargeAV_outcorr.gpl'
corry = 'onemp_y/onempy_finephys_bph_-m432.gpl'
corrz = 'onemp_z/onempz_finephys_bph_-m432.gpl'

#corrpathx  = os.path.join(basedir, ensemble, corrx)
corrpathy  = os.path.join(basedir, ensemble, corry)
corrpathz  = os.path.join(basedir, ensemble, corrz)

#outfile = mass+'_Rho_iAV_Q0.202_chargeAV_outcorr.gpl'
outfile = 'onemp_finephys_bph_-m432.gpl'
outpath = os.path.join(basedir, ensemble, outfile)

#tags = ['onemm.HH', 'onemm.Hh', 'onemm.HR', 'onemm.Hr', 'onemm.RR', 'onemm.RH', 'onemm.Rr', 'onemm.Rh', 'onemm.hH', 'onemm.hR', 'onemm.hr', 'onemm.hh', 'onemm.rR', 'onemm.rH', 'onemm.rh', 'onemm.rr']
tags=['onemp.ll', 'onemp.gl', 'onemp.lg', 'onemp.gg']

#corrx = gv.dataset.Dataset(corrpathx) #, grep='onempx')
corry = gv.dataset.Dataset(corrpathy) #, grep='onempy')
corrz = gv.dataset.Dataset(corrpathz) #, grep='onempz')

corrav = []

#xkeys=corrx.keys() 
ykeys=list(corry.keys()) 
zkeys=list(corrz.keys()) 

for i in range(len(corry)):
  tmp=np.average((np.array([corry[ykeys[i]], corrz[zkeys[i]]])), axis=0)
  #tmp=np.average((np.array([corrx[xkeys[i]], corry[ykeys[i]], corrz[zkeys[i]]])), axis=0)
  print(tmp.shape)
  corrav.append(tmp)

corr_av  = gv.dataset.avg_data(corrav[0])	
#corrx_av = gv.dataset.avg_data(corrx)		# just one pol
corry_av = gv.dataset.avg_data(corry)		# just one pol
corrz_av = gv.dataset.avg_data(corrz)		# just one pol

#tst = corr_av/(corrx_av+corry_av+corrz_av)*3

diff = gv.sdev(corr_av)/gv.sdev(corry_av[ykeys[0]])
print(diff)
print(type(corrav[0][0]))

ct = 0
with open(outpath, 'w') as f:
	for cfgs in corrav:			
		for cfg in cfgs:
			f.write(tags[ct]+'  ')			
			cfg.tofile(f, sep=" ", format="%.17g")
			f.write('\n')
		ct = ct + 1

