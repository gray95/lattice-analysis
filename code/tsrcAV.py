# averaging over multiple time sources

import gvar as gv
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

# l3248f211b580m002426m06730m8447
# l3264f211b600m00507m0507m628a-coul-v5
# l3296f211b630m0074m037m440-coul-v5
# l4864f211b600m001907m05252m6382 - coarse phys
# l6496f211b630m0012m0363m432 - fine phys
ensemble = 'l6496f211b630m0012m0363m432/l6496f211b630m0012m0363m432h-coul'
basedir = '../data/hybrid'

outfile = 'tav_onempx_fine_m432.gpl' 

t0 = [y for y in range(0,96,6)]

corr = ['t'+str(i)+'_onemp_fine_m432.gpl' for i in t0]
print(corr)

#tags=['onemmy.HH', 'onemmy.Hh', 'onemmy.HR', 'onemmy.Hr', 'onemmy.RR', 'onemmy.RH', 'onemmy.Rr', 'onemmy.Rh', 'onemmy.hH', 'onemmy.hR', 'onemmy.hr', 'onemmy.hh', 'onemmy.rR', 'onemmy.rH', 'onemmy.rh', 'onemmy.rr']
#tags=['onempx.ll', 'onempy.ll', 'onempz.ll', 'onempx.gl', 'onempy.gl', 'onempz.gl', 'onempx.lg', 'onempy.lg', 'onempz.lg', 'onempx.gg', 'onempy.gg', 'onempz.gg']
tags=['onempx.ll', 'onempx.gl', 'onempx.lg', 'onempx.gg']

corrpath = []
C        = []
for corr in corr:
	C.append(gv.dataset.Dataset(os.path.join(basedir, ensemble, corr)))

corrav = []
keys = list(C[0].keys())
print(type(keys))
#sys.exit(0)
for key in keys:
  tmp=np.average( (np.array( [C[j][key] for j in range(len(C))] ) ) , axis=0) 
  print(tmp.shape)
  corrav.append(tmp)


corr_av_ll=gv.dataset.avg_data(corrav[0])

corrt0_av = gv.dataset.avg_data(C[0])

diff = gv.sdev(corr_av_ll)/gv.sdev(corrt0_av[keys[0]])

print(diff)
outpath = os.path.join(basedir, ensemble, outfile)



ct = 0
with open(outpath, 'w') as f:
	for cfgs in corrav:			
		for cfg in cfgs:
			f.write(tags[ct]+'  ')
			cfg.tofile(f, sep=" ", format="%.10g")
			f.write('\n')
		ct = ct + 1

