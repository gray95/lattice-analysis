import numpy as np
import matplotlib.pyplot as plt

y        = dict()
amus     = dict()
amus_err = dict()
sy       = dict()

lbmw = "BMW 2020"
lrbc = "RBC/UKQCD 2018"
letm = "ETM 2020"
lus  = "this work"

y[lbmw]        = 3
amus[lbmw]     = 0.0086
amus_err[lbmw] = 0.0059
sy[lbmw]  =  "ko"

y[lrbc]        = 4
amus[lrbc]     = 0.0149
amus_err[lrbc] = 0.0032
sy[lrbc]       =  "ko"

y[letm]        = 5
amus[letm]     = 0.0053
amus_err[letm] = 0.0033
sy[letm]       =  "ko"

#-2.63(41)e-12
y[lus]        = 6
amus[lus]     = 0.0263
amus_err[lus] = 0.0041
sy[lus]       =  "k*"

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

plt.figure(figsize=(5,6))

for kkk in y :
   plt.errorbar(x=amus[kkk], xerr=amus_err[kkk] , y=y[kkk], fmt=sy[kkk])
   if kkk==lus:
      plt.text(x=amus[kkk]-0.5*amus_err[kkk], y=0.95*y[kkk], s=kkk)
   else:
      plt.text(x=amus[kkk]-0.5*amus_err[kkk], y=1.05*y[kkk], s=kkk)



plt.xlabel(r'$- \delta  a_{\mu}^{s} \times 10^{10}$', fontsize=18)

plt.gca().axes.get_yaxis().set_visible(False)
plt.tight_layout()
plt.savefig("../figures/amusdiff_summary.png", dpi=500)

plt.show()

