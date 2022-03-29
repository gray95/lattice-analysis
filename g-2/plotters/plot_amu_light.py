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

#-4.0(1.4)e-11
tagUS='this work (preliminary)'
y[tagUS]        = 4
amu_light[tagUS]     = -0.40
amu_light_err[tagUS] =  0.14
sy[tagUS]       =  "ko"


for kkk in y :
   plt.errorbar(x=amu_light[kkk], xerr=amu_light_err[kkk] , y=y[kkk], fmt=sy[kkk] , label = kkk)

plt.plot([0,0],[0,100],"k--", alpha=0.5)

plt.title(r"(QCD+QED-QCD)  $\delta a_{\mu}^{light}$")
plt.legend(loc="lower right")


plt.xlabel(r'$\delta  a_{\mu}^{light} \times 10^{10}$')

plt.xlim(-2, 2)
plt.ylim(1, 8)

plt.gca().axes.get_yaxis().set_visible(False)

plt.savefig("../figures/amu_light_summary.png")

plt.show()

