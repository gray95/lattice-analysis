import numpy as np
import matplotlib.pyplot as plt


y        = dict()
amu_light     = dict()
amu_light_err = dict()
sy        = dict()

tagA='BMW Nature'
y[tagA]        = 7
amu_light[tagA]     = -1.28
amu_light_err[tagA] = 0.59
sy[tagA]  =  "ro"

tagB="RBC/UKQCD 1801.07224"
y[tagB]        = 6
amu_light[tagB]     = 0.2
amu_light_err[tagB] = 0.2
sy[tagB]       =  "go"

tagC="ETMC 1901.10462"
y[tagC]        = 5
amu_light[tagC]     = 1.1
amu_light_err[tagC] = 1.0
sy[tagC]       =  "bo"

tagUS='coarse-1.5fm-scipy'
y[tagUS]        = 4
amu_light[tagUS]     = -0.136
amu_light_err[tagUS] =  0.06
sy[tagUS]       =  "ko"
tagUS1='coarse-2fm-gvar'
y[tagUS1]        = 3
amu_light[tagUS1]     = -0.14
amu_light_err[tagUS1] =  0.28
sy[tagUS1]       =  "kx"
tagUS1='coarse-1.5fm-gvar'
y[tagUS1]        = 2
amu_light[tagUS1]     = -0.21
amu_light_err[tagUS1] =  0.24
sy[tagUS1]       =  "k+"





for kkk in y :
   plt.errorbar(x=amu_light[kkk], xerr=amu_light_err[kkk] , y=y[kkk], fmt=sy[kkk] , label = kkk)

plt.plot([0,0],[0,100],"k--", alpha=0.5)

plt.title(r"(QCD+QED-QCD)  $\delta a_{\mu}^{light}$")
plt.legend(loc="lower right")


plt.xlabel(r'$\delta  a_{\mu}^{light} \times 10^{10}$')

plt.xlim(-2, 2)
plt.ylim(1, 8)

plt.gca().axes.get_yaxis().set_visible(False)

#plt.savefig("amu_light_summary.png")

plt.show()

