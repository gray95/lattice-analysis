
import matplotlib.pyplot as plt
import numpy as np
import gvar as gv
import sys

#Y = gv.load('../fits/tstar_fine.p')
Y = gv.load('../fits/tstar_fine_3ml.p')

original = Y['corr']
new      = Y['newcorr']
new_t    = {}
ratio = {}

otags = list(original.keys())
ntags = list(new.keys())

## fold reconstructed corr to make them periodic 
for K in ntags:
    new_t[K]= 0.5*( new[K] + np.flip(new[K]) )

for (t1,t2) in zip(otags,ntags):
    ratio[t2] = new_t[t2]/original[t1]       #ratio[t2] = new[t2]/original[t1]
T = len(original[otags[0]])
print(ratio)

#############################################################################

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

fig, ax1 = plt.subplots()

#ax1.errorbar(range(T), [y.mean for y in ratio['nocharge']], xerr=0, yerr=[y.sdev for y in ratio['nocharge']], fmt='h', label='corr', color='blue')
ax1.errorbar(range(T), [y.mean for y in new_t['nocharge']], xerr=0, yerr=[y.sdev for y in new_t['nocharge']], fmt='h', label='reconstruction with folding', color='red')
#ax1.errorbar(range(T), [y.mean for y in new['nocharge']], xerr=0, yerr=[y.sdev for y in new['nocharge']], fmt='h', label='reconstruction', color='red')
ax1.errorbar(range(T), [y.mean for y in original[otags[0]]], xerr=0, yerr=[y.sdev for y in original[otags[0]]], fmt='h', label='data', color='blue')

ax1.axvline(x=23, linestyle='--', label='$t^*$')
#ax1.set_xlim(left=20,right=48)
plt.legend(frameon=False,fontsize=8,loc='lower right')
plt.yscale('log')

plt.tight_layout()
#plt.savefig('../figures/tstar_recon_test.png', dpi=500, bbox_inches="tight")

plt.show()
