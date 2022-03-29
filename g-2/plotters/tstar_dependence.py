import gvar as gv
import sys
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
import g2tools as g2
import matplotlib.pyplot as plt
import numpy as np

## ensemble
a_str = 'vc'

a = (w0/w0overa[a_str])/hbarc           # in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV


## load in corrs
C = gv.load('../fits/vc_rhocorr.p')

mq = ["3ml", "5ml", "7ml"]
keys = ["VTvc3mlq0","VTvc5mlq0","VTvc7mlq0"]
data = C['3ml']
T = int(len(data['nocharge']))
tcut = T
tstar = 13
Y = {}

## calculate amu for different values of tstar
for m,M in zip(mq,keys):
    print("mq = %s"%m)
    Y[m] = []
    original_data = C['org'][M]
    for t in range(tstar, T):
        ddata = original_data[:t]
        ddata = np.append(ddata, C[m]['nocharge'][t:])
        print(len(ddata))
        print("t* = %d"%t)
        vpol = g2.fourier_vacpol(ddata, Z=ZV, ainv=1/a, periodic=False)
        unchargedamud = g2.a_mu(vpol,1/3.)
        unchargedamuu = g2.a_mu(vpol,2/3.)

        amu = unchargedamud + unchargedamuu
        Y[m].append(amu)

print(Y)
print(T)

plt.errorbar(range(tstar,T), [y.mean for y in Y['3ml']], yerr=[y.sdev for y in Y['3ml']], fmt='h', color='red', label='3ml')
#plt.errorbar(range(1,T), [y.mean for y in Y['5ml']], yerr=[y.sdev for y in Y['5ml']], fmt='h', color='blue', label='5ml')
#plt.errorbar(range(1,T), [y.mean for y in Y['7ml']], yerr=[y.sdev for y in Y['7ml']], fmt='h', color='green', label='7ml')
plt.title('very coarse ensemble')
plt.xlabel('tstar')
plt.ylim(bottom=4.4e-08, top=6.5e-08)
plt.xlim(left=tstar, right=36)
plt.ylabel('amu')
#plt.savefig('../figures/tcut_dependence.png')
plt.show()


