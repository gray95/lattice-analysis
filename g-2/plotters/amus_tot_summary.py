import numpy as np
import gvar as gv
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from commonweal import hpqcd_amus, bmw_amus

y        = dict()
amus     = dict()
amus_err = dict()
sy       = dict()

lbmw = 'BMW 2020'
lhpqcd = 'HPQCD 2014'
lrbc = 'RBC/UKQCD 2018'
letm = 'ETM 2020'
lus  = 'this work (blinded)'
print(hpqcd_amus)
print(bmw_amus)
y[lbmw]        = 5
amus[lbmw]     = bmw_amus.mean
amus_err[lbmw] = bmw_amus.sdev
sy[lbmw]  =  "ko"

y[lhpqcd]        = 2
amus[lhpqcd]     = hpqcd_amus.mean
amus_err[lhpqcd] = hpqcd_amus.sdev
sy[lhpqcd]       =  "ko"

y[lrbc]        = 4
amus[lrbc]     = 53.2e-10
amus_err[lrbc] = 0.4e-10
sy[lrbc]       =  "ko"

y[letm]        = 3
amus[letm]     = 53.1e-10
amus_err[letm] = 2.5e-10
sy[letm]       =  "ko"

y[lus]        = 6
amus[lus]     = 54.38e-10
amus_err[lus] = 0.80e-10
sy[lus]       =  "k*"

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

plt.figure(figsize=(4,6))

for kkk in y :
   plt.errorbar(x=amus[kkk], xerr=amus_err[kkk] , y=y[kkk], fmt=sy[kkk])
   if kkk==lus:
       plt.text(x=amus[kkk]-2.0*amus_err[kkk], y=0.97*y[kkk], s=kkk)
   else:
       plt.text(x=amus[kkk]-0.7*amus_err[kkk], y=1.03*y[kkk], s=kkk)
       

plt.xlabel(r'$a_{\mu}^{(s)}$', fontsize=20, labelpad=10)

plt.gca().axes.get_yaxis().set_visible(False)
plt.tight_layout()
plt.savefig("../figures/amustot_summary.png", dpi=500)

plt.show()

