import numpy as np
import matplotlib.pyplot as plt

y        = dict()
amus     = dict()
amus_err = dict()
sy       = dict()

lbmw = "BMW 2020"
lrbc = "RBC/UKQCD 2018"
letm = "ETM 2020"
lus  = "this work (blinded)"

y[lbmw]        = 3
amus[lbmw]     = 0.86
amus_err[lbmw] = 0.59
sy[lbmw]       = "ko"

y[lrbc]        = 4
amus[lrbc]     = 1.49
amus_err[lrbc] = 0.32
sy[lrbc]       = "ko"

y[letm]        = 5
amus[letm]     = 0.53
amus_err[letm] = 0.33
sy[letm]       = "ko"

y[lus]        = 6
amus[lus]     = 0.92
amus_err[lus] = 0.81
sy[lus]       = "k*"

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

plt.figure(figsize=(4,6))

for kkk in y :
   plt.errorbar(x=amus[kkk], xerr=amus_err[kkk] , y=y[kkk], fmt=sy[kkk])
   plt.text(x=amus[kkk]-1.0*amus_err[kkk], y=1.03*y[kkk], s=kkk)



plt.xlabel(r'$- \delta  a_{\mu}^{s} \times 10^{12}$', fontsize=18)
plt.xlim(left=-0.005)
plt.gca().axes.get_yaxis().set_visible(False)
plt.tight_layout()
plt.savefig("../figures/amusdiff_summary.png", dpi=500)

plt.show()

