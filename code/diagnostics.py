# VARIOUS SMALL METHODS TO GET A FEEL FOR THE DATA
# compute how closely two different correlators are correlated
# compute the fractional error in the correlator

import gvar as gv
import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

bp = '../data/hybrid'
conf= ''
#l3296f211b630m0074m037m440-coul-v5
#l3264f211b600m00507m0507m628a-coul-v5
#l3248f211b580m002426m06730m8447
#l4864f211b600m001907m05252m6382
gpls = ['l3248f211b580m002426m06730m8447/onemp_vc_abij_m863.gpl', ] 

files = [os.path.join(bp, conf, gpl) for gpl in gpls]

corrs = [gv.dataset.Dataset(f) for f in files]

corrs = [gv.dataset.avg_data(c) for c in corrs]

corr1 = corrs[0]['onemp.ll']
#corr2 = corrs[1]['onemm.HH']
#corr3 = corrs[1]['onemm.RR']
#corr4 = corrs[1]['onemp.gg']
#corr5 = corrs[4]['onemm.RR']

y1 = [abs(x.sdev/x.mean) for x in corr1]
#y2 = [abs(x.sdev/x.mean) for x in corr2]
#y3 = [abs(x.sdev/x.mean) for x in corr3]
#y4 = [abs(x.sdev/x.mean) for x in corr4]
#y5 = [abs(x.sdev/x.mean) for x in corr5]

print(corr1)
print(y1)

sys.exit(0)


## fit data to this function
def f( t, A, B):
		return A*np.exp(B*t)

T = 9	# stop plotting here
TT = 7  # stop fitting here

## do the fit
#start = (0.1, 0.5)		# parameter search starts here
#popt, pcov = curve_fit( f, range(TT), y1[:TT], p0=None)
#print(popt)
## LEPAGE arg
#M_ps = gv.gvar('1.366839(72)')
#M_h  = gv.gvar('1.934(80)')
#M_h  = gv.gvar('1.934(80)')
#print(M_ps.mean - 0.5*M_h.mean)
#lepage = [ popt[0]*(gv.exp((M_h - M_ps)*t)).mean for t in range(TT) ]
#fit_exp = [ popt[0]*(np.exp(popt[1]*t)) for t in range(TT) ]

sys.exit(0)

fig, ax = plt.subplots()

#plt.ylim(top=1)
plt.yscale('log')
plt.xlabel('t/a')
plt.ylabel(r'$\frac{\sigma_G}{\overline{G}}', rotation=0, labelpad = 20)
ax.yaxis.set_label_coords(-.1, 0.43)
#plt.plot([(y1[i]/y3[i]) for i in range(12)], 'ro', label='32 over 10 ape')
plt.plot(y1[:T], 'ro', label=r'$1^{-+}$ hybrid')
#plt.plot(lepage,  'r--', label='$1^{-+}$ lepage argument')
#plt.plot(y2[:T], 'gx', label=r'$1^{--}$ hybrid')
#plt.plot(y3[:T], 'b+', label=r'$\bar{\psi}\gamma_i\psi$', markersize=10)
#plt.plot(fit_exp, 'rx', label='fit to exp')
#plt.plot(y5[:15], 'bo', label=r'$1^{--} \bar{\psi}\gamma_i\psi$')
plt.title('error/mean of correlator', rotation=0)
#plt.title(conf)
plt.text(5, 0.0003, 'FINE ENSEMBLE', fontsize=15)
plt.legend(fontsize=12)
#plt.savefig('../figures/fractional_error_corr_fineh.png', dpi=500, bbox_inches="tight")
plt.show()
