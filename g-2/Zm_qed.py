
## fit equation 15 to data in Table VIII of the charmonium paper [arXiv:2005.01845]
## extrapolate to coarser lattice - vc specifically.
## need this for quark mass renorm in ETM scheme with QED

# goal: fit RQED[Zm] to f((mu*a)**2)
import gvar as gv
import numpy as np
from gvar import log
from commonweal import w0, hbarc
import sys
import lsqfit as lsq
# start with data

w0overa = {'S5':gv.gvar('1.4029(9)'), 'S10':gv.gvar('1.9330(20)'), 'S12':gv.gvar('2.8941(48)')}

a = [ w0/(hbarc*w0overa[x]) for x in w0overa ] 
print("lattice spacings in 1/GeV: %s"%a)

Rqed_Zm_SMOM = gv.BufferDict()

Rqed_Zm_SMOM['S5'] = {2:gv.gvar('1.001200(83)'), 2.5:gv.gvar('1.000827(31)'), 3:gv.gvar('1.000540(15)')}
Rqed_Zm_SMOM['S10'] = {2:gv.gvar('1.001516(35)'), 3:gv.gvar('1.000851(11)'), 4:gv.gvar('1.0005001(21)')}
Rqed_Zm_SMOM['S12'] = {2:gv.gvar('1.001853(83)'), 3:gv.gvar('1.001308(18)'), 4:gv.gvar('1.0009331(34)')}

corrs5 = {}     # none listed so corrs from set 2
corrs10 = {}
corrs12 = {}
## need to put in correlations between different mu values on same lattice
# gv.correlate()


corrs5[2,2] = corrs5[2.5,2.5] = corrs5[3,3] = 1
corrs5[2,2.5] = corrs5[2.5,2] = 0.41
corrs5[2,3] = corrs5[3,2] = 0.12
corrs5[2.5,3] = corrs5[3,2.5] = 0.45
corrs10[2,2] = corrs10[3,3] = corrs10[4,4] = 1
corrs10[2,3] = corrs10[3,2] = 0.41
corrs10[2,4] = corrs10[4,2] = 0.52
corrs10[3,4] = corrs10[4,3] = 0.42
corrs12[2,2] = corrs12[3,3] = corrs12[4,4] = 1
corrs12[2,3] = corrs12[3,2] = 0.26
corrs12[2,4] = corrs12[4,2] = 0.19
corrs12[3,4] = corrs12[4,3] = 0.26


Rqed_Zm_SMOM['S5'] = gv.correlate(Rqed_Zm_SMOM['S5'], corrs5, upper=True, verify=True)
Rqed_Zm_SMOM['S10'] = gv.correlate(Rqed_Zm_SMOM['S10'], corrs10, upper=True, verify=True)
Rqed_Zm_SMOM['S12'] = gv.correlate(Rqed_Zm_SMOM['S12'], corrs12, upper=True, verify=True)


Zm_MS_SMOM = {2:0.999872, 2.5:0.999873, 3:0.999873, 4:0.999873}
rqed_3gev = {2:0.999372, 2.5:0.999718, 3:1, 4:1.000446}

Rqed_Zm_MS = []

for S in Rqed_Zm_SMOM:
    for mu in Rqed_Zm_SMOM[S].keys():
        Rqed_Zm_MS.append(Rqed_Zm_SMOM[S][mu]*Zm_MS_SMOM[mu])#*rqed_3gev[mu])

aqed = 1/134   # at charm mass
Q     = 2/3
#aqcd  = {2:0.5 , 2.5:0.45 , 3:0.4,  4:0.35}
aqcd  = {2:gv.gvar('0.3030(54)') , 2.5:gv.gvar('0.2741(43)') , 3:gv.gvar('0.2545(37)'),  4:gv.gvar('0.2291(29)')}
K     = aqed*Q**2     
EXT=False

def fcn(p):
    Zm = []
    A  = p['a']
    sets = ['S5','S10','S12']
    C = p['C'] 
    for k in sets:
        MU = Rqed_Zm_SMOM[k].keys()
        if EXT:
            a = A
            MU = p['MU']
        elif k=='S5':
            a = A[0]
        elif k=='S10':
            a = A[1]
        elif k=='S12':
            a = A[2]
        for mu in MU:
            D = U = 0
            u = mu
            u2 = (1/u)**2
            di = (u*a)**2
            #di = 2*log(u*a)
            for j in range(len(p['di'])):
                D += p['di'][j]*di**(j+1)
            for n in range(len(p['mu'])):
                U += p['mu'][n]*u2**(n+1)
            Y = C * (1 + K * ( D + aqcd[mu].mean*U ) )
            Zm.append( Y )
    return Zm

#for Nd in range(1,5):
#for Nu in range(1,5):
Nd = 3
Nu = 1
prior = gv.BufferDict()
prior['C'] = gv.gvar('1.0(1)')
prior['di'] = [gv.gvar('0(1)') for n in range(Nd)]
prior['mu'] = [gv.gvar('0(1)') for n in range(Nu)]
prior['a'] = a
fit = lsq.nonlinear_fit(data=Rqed_Zm_MS, prior=prior, fcn=fcn)
if fit.chi2/fit.dof < 1:
    print("%d a**2 and %d 1/u**2"%(Nd,Nu))
    print(fit.chi2/fit.dof)
print(fit)
#print(prior)
print(fit.p)
## extrapolate to very coarse
EXT=True
X = fit.p
#X['C'] = fit.p['C']
#X['mu'] = fit.p['mu']
#X['di'] = fit.p['di']
X['a']  = gv.gvar('0.7677(40)')
X['MU'] = [2]
print(X)
#0.6195(33)
# vc - 0.7677(40)
print("\nEXTRAPOLATED VALUES AT VARIOUS scales")
ext = fcn(X)
print(ext)
