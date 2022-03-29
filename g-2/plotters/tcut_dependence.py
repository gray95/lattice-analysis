
import gvar as gv
import sys
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
import g2tools as g2
import matplotlib.pyplot as plt
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

data = C['3ml']
T = int(len(data['nocharge']))
Y = {}
## calculate amu for different values of tcut

for m in mq:
    data = C[m]
    Y[m] = []
    for t in range(1,T):
        print("tcut = %d"%t)
        vpol = g2.fourier_vacpol(data['nocharge'][:t], Z=ZV, ainv=1/a, periodic=False)
        unchargedamud = g2.a_mu(vpol,1/3.)
        unchargedamuu = g2.a_mu(vpol,2/3.)

        vpol = g2.fourier_vacpol(data['down-charge'][:t], Z=ZVqedd, ainv=1/a, periodic=False)
        chargedamud = g2.a_mu(vpol,1/3.)

        vpol = g2.fourier_vacpol(data['up-charge'][:t], Z=ZVqedu, ainv=1/a, periodic=False)
        chargedamuu = g2.a_mu(vpol,2/3.)

        amu = unchargedamud + unchargedamuu
        Y[m].append(amu)

## plot amu v tcut

plt.errorbar(range(1,T), [y.mean for y in Y['3ml']], yerr=[y.sdev for y in Y['3ml']], fmt='h', color='red', label='3ml')
plt.errorbar(range(1,T), [y.mean for y in Y['5ml']], yerr=[y.sdev for y in Y['5ml']], fmt='h', color='blue', label='5ml')
plt.errorbar(range(1,T), [y.mean for y in Y['7ml']], yerr=[y.sdev for y in Y['7ml']], fmt='h', color='green', label='7ml')
plt.title('very coarse ensemble')
plt.xlabel('tcut')
plt.ylabel('amu')
plt.savefig('../figures/tcut_dependence.png')
plt.show()
