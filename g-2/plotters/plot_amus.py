import numpy as np
import matplotlib.pyplot as plt


y        = dict()
amus     = dict()
amus_err = dict()
sy        = dict()

y['BMW']        = 3
amus['BMW']     = 0.0086
amus_err['BMW'] = 0.0050
sy['BMW']  =  "ro"

y['RBC/UKQCD']        = 4
amus['RBC/UKQCD']     = 0.0149
amus_err['RBC/UKQCD'] = 0.0032
sy['RBC/UKQCD']       =  "go"

y['ETMC']        = 5
amus['ETMC']     = 0.0053
amus_err['ETMC'] = 0.0030
sy['ETMC']       =  "bo"

#-2.63(41)e-12
y['This work (preliminary)']        = 6
amus['This work (preliminary)']     = 0.0263 
amus_err['This work (preliminary)'] = 0.0041
sy['This work (preliminary)']       =  "ko"



for kkk in y :
   plt.errorbar(x=amus[kkk], xerr=amus_err[kkk] , y=y[kkk], fmt=sy[kkk] , label = kkk)



plt.title(r"(QCD+QED-QCD)  $\delta a_{\mu}^{s}$")
plt.legend(loc='lower right')


plt.xlabel(r'$- \delta  a_{\mu}^{s} \times 10^{10}$')

plt.xlim(0, 0.04)

plt.gca().axes.get_yaxis().set_visible(False)

plt.savefig("../figures/amu_s_summary.png")

plt.show()

