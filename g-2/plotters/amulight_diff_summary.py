import numpy as np
import matplotlib.pyplot as plt

y        = dict()
amu_light     = dict()
amu_light_err = dict()
sy        = dict()

lbmw = "BMW 2020"
lrbc = "RBC/UKQCD 2018"
letm = "ETM 2020"
lus  = "this work"

y[lbmw]        = 3
amu_light[lbmw]     = -1.23
amu_light_err[lbmw] = 0.50
sy[lbmw]  =  "ko"

y[lrbc]        = 1
amu_light[lrbc]     = 5.9
amu_light_err[lrbc] = 5.9
sy[lrbc]       =  "ko"

y[letm]        = 2
amu_light[letm]     = 1.1
amu_light_err[letm] = 1.0
sy[letm]       =  "ko"

y[lus]        = 4
amu_light[lus]     = -0.4
amu_light_err[lus] =  0.4
sy[lus]       =  "k*"


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

plt.figure(figsize=(5,6))

for kkk in y :
   plt.errorbar(x=amu_light[kkk], xerr=amu_light_err[kkk] , y=y[kkk], fmt=sy[kkk])
   if kkk==lus:
      plt.text(x=amu_light[kkk]-0.5*amu_light_err[kkk], y=0.95*y[kkk], s=kkk)
   else:
      plt.text(x=amu_light[kkk]-0.5*amu_light_err[kkk], y=1.05*y[kkk], s=kkk)

plt.axvline(x=0, linestyle='--', alpha=0.5)

plt.xlabel(r'$\delta  a_{\mu}^{(l)} \times 10^{10}$', fontsize=20)

plt.gca().axes.get_yaxis().set_visible(False)

plt.tight_layout()
plt.savefig("../figures/amu_light_summary.png", dpi=500)

plt.show()

