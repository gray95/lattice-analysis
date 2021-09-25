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
conf= 'l3296f211b630m0074m037m440-coul-v5'
#l3296f211b630m0074m037m440-coul-v5
#l3264f211b600m00507m0507m628a-coul-v5
#l3248f211b580m002426m06730m8447
#l4864f211b600m001907m05252m6382
gpls = ['onemp_fine_m450.gpl'] 

files = [os.path.join(bp, conf, gpl) for gpl in gpls]

corrs = [gv.dataset.Dataset(f) for f in files]

corrs = [gv.dataset.avg_data(c) for c in corrs]

corr1 = corrs[0]['onemp.ll']
#corr2 = corrs[1]['onemmy.hh']
#corr3 = corrs[1]['onemp.ll']
#corr4 = corrs[1]['onemp.gg']
#corr5 = corrs[4]['onemm.RR']

y1 = [abs(x.sdev/x.mean) for x in corr1]
#y2 = [abs(x.sdev/x.mean) for x in corr2]
#y3 = [abs(x.sdev/x.mean) for x in corr3]
#y4 = [abs(x.sdev/x.mean) for x in corr4]
#y5 = [abs(x.sdev/x.mean) for x in corr5]

## fit data to this function
def f( t, A, B):
		return A*np.exp(B*t)

T = 7	# stop fitting/plotting here

## do the fit
#start = (0.1, 0.5)		# parameter search starts here
popt, pcov = curve_fit( f, range(T), y1[:T], p0=None)
print(popt)
print(pcov)

## LEPAGE arg
M_ps = gv.gvar('1.366839(72)')
#M_h  = gv.gvar('1.934(80)')
M_h  = gv.gvar('2.000(80)')
print(M_ps.mean - 0.5*M_h.mean)
lepage = [ 0.0038*(gv.exp((2*M_ps - M_h)*0.5*t)).mean for t in range(T) ]
fit_exp = [ popt[0]*(np.exp(popt[1]*t)) for t in range(T) ]

plt.yscale('log')
plt.xlabel('t/a')
plt.ylabel('fractional\n error', rotation=0)
#plt.plot([(y1[i]/y3[i]) for i in range(12)], 'ro', label='32 over 10 ape')
plt.plot(y1[:T], 'ro', label=r'$1^{-+}$ hybrid data')
plt.plot(lepage,  'b--', label='lepage argument')
plt.plot(fit_exp, 'gx', label='fit to exp')
#plt.plot(y2[:12], 'go', label='16 src')
#plt.plot(y4[:15], 'g+', label='gg-10ape')
#plt.plot(y5[:15], 'bo', label=r'$1^{--} \bar{\psi}\gamma_i\psi$')
plt.title('error/mean of correlator')
#plt.title(conf)
plt.text(3, 0.004, 'FINE ENSEMBLE', fontsize=10)
plt.legend()
plt.savefig('../figures/fractional_error_corr_lepage.png', dpi=500, bbox_inches="tight")
plt.show()
