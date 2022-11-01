
# checking that mu and md are done through \Delta M^2

import sys
import matplotlib
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from commonweal import w0overa, w0, hbarc
from gvar import log

## IDEA: fit M^2(Q=1/3) and M^2(Q=2/3) to a func of m_q and a^2
a_str = ['vc', 'c', 'f']

PS_sqmass_n = []
PS_sqmass_u = []
PS_sqmass_d = []
for a in a_str:
    if a=='vc':
        ll = 'ps_vcoarse'
    elif a=='c':
        ll = 'ps_coarse'
    elif a=='f':
        ll = 'ps_fine'

    obs = gv.load('./pseudoscalar/fits/'+ll+'.p')
    overa = (w0overa[a]/w0)*hbarc

    PS_sqmass_n.extend( [ (y*overa)**2 for y in obs['E:n'] ] )
    PS_sqmass_d.extend( [ (y*overa)**2 for y in obs['E:d'] ] )
    PS_sqmass_u.extend( [ (y*overa)**2 for y in obs['E:u'] ] )

    del PS_sqmass_n[-1]
    del PS_sqmass_u[-1]
    del PS_sqmass_d[-1]

    break


a = [ w0/(hbarc*w0overa[s]) for s in a_str ]
xdata = {'a':a}
ydata = {'up':PS_sqmass_u, 'down':PS_sqmass_d }
print(ydata)
##########################################################################

def fcn_PSmass2(p):
    out = {'up':[], 'down':[]}
    A = p['a']
    M = [3, 5, 7]
    u0, u1, u2, u3, u4 = p['uC']
    d0, d1, d2, d3, d4 = p['dC']
#    for a in A:
    for m in M:
        out['up'].append(   u0 * (1 + u1*m + u2*m*log(m) + u4*m**2*log(m)) )
        out['down'].append( d0 * (1 + d1*m + d2*m*log(m) + d4*m**2*log(m)) )
    return out

def make_prior_PSmass2(x):
    prior = gv.BufferDict()
    L = 10
    N = 4
    prior['uC'] = gv.gvar(['0(1)e-3']+[gv.gvar(0,L) for n in range(N)])
    prior['dC'] = gv.gvar(['0(1)e-3']+[gv.gvar(0,L) for n in range(N)])
    prior['a'] = x['a']
    return prior

##########################################################################

print("Fitting squared PS mass")
print("------------------------------------------------------")
uprior = make_prior_PSmass2(xdata)
print(uprior)
ufit = lsq.nonlinear_fit(data=ydata,prior=uprior,fcn=fcn_PSmass2)
print(ufit)

