
# This script calculates the additive factors needed for a Dashen like EM renorm scheme
#  -damu/dmq * mq * delta_u/d(mq)
# first term from a fit to amu data with a spline
# third term from chiral fits to PS data

import sys
import matplotlib
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from commonweal import w0overa, w0, hbarc
from gvar import log

vc_vt = 'out/vc_no_vector_fit.p'
c_vt  = 'out/c_no_vector_fit.p'
f_vt  = 'out/f_no_vector_fit.p'

vc_ps = 'pseudoscalar/fits/ps_vcoarse.p' 
c_ps  = 'pseudoscalar/fits/ps_coarse.p' 
f_ps  = 'pseudoscalar/fits/ps_fine.p' 

a_str = 'vc'

if a_str=='vc':
    vt = gv.load(vc_vt)
    ps = gv.load(vc_ps)
elif a_str=='c':
    vt = gv.load(c_vt)
    ps = gv.load(c_ps)
elif a_str=='f':
    vt = gv.load(f_vt)
    ps = gv.load(f_ps)
else:
    print("something went awry")
###################################################################################
####### fitting functions + priors ################
def PSfcn(x,p):
    res = []
    for pt in x:
        res.append( p[0]*(1 + p[1]*pt + p[2]*pt**2 ) )
    return res

def make_prior(x):
    prior = gv.BufferDict()
    prior['uC'] = gv.gvar(['4(2)e-8','0.1(10)','0.1(10)','0.1(10)','0.1(10)'])
    prior['dC'] = gv.gvar(['1(1)e-8','0.1(10)','0.1(10)','0.1(10)','0.1(10)'])
    prior['mq']= x['mq']
    return prior

def amufcn(p):
    out = {'up':[], 'down':[]}
    M = p['mq']
    u0, u1, u2, u3, u4 = p['uC']
    d0, d1, d2, d3, d4 = p['dC']
    #u2 = d2 = 0
    u3 = d3 = 0
    u4 = d4 = 0
    for m in M:
        out['up'].append( u0 * (1 + u1*m + u2*m*log(m) + u3*m**1.5 + u4*m**2 ) )  
        out['down'].append( d0 * (1 + d1*m + d2*m*log(m) + d3*m**1.5 + d4*m**2 ) ) 
    return out
###############################################################################
# PS fits
am = [3, 5, 7]
overa = (w0overa[a_str]/w0)*hbarc
obs = ps
PS_sqmass_n = [ y**2 for y in obs['E:n'] ]
PS_sqmass_d = [ y**2 for y in obs['E:d'] ]
PS_sqmass_u = [ y**2 for y in obs['E:u'] ]
del PS_sqmass_n[-1]
del PS_sqmass_u[-1]
del PS_sqmass_d[-1]

del_uu = [ (A**2-C**2) for (A,C) in zip(obs['E:u'],obs['E:n']) ]
del_dd = [ (B**2-C**2) for (B,C) in zip(obs['E:d'],obs['E:n']) ]
del del_uu[-1]
del del_dd[-1]

xdata = am #np.array([ x for x in am ])
ydata = del_uu
y_data = [ y.mean for y in del_uu ]
y_err  = [ y.sdev for y in del_uu ]

L = 500
print("Fitting UP squared diff")
print("------------------------------------------------------")
uprior = gv.gvar([gv.gvar('0(1)e-04'), gv.gvar(0,L), gv.gvar('0(1)')])
ufit = lsq.nonlinear_fit(data=(xdata, ydata),prior=uprior,fcn=PSfcn)
print(ufit)
print("NOISY FIT")
noisyufit = lsq.nonlinear_fit(data=(xdata, ydata),prior=uprior,fcn=PSfcn, noise=True)
print("noisy reduced chi2 is %s"%(noisyufit.chi2/noisyufit.dof))

print("Fitting DOWN squared diff")
print("------------------------------------------------------")
ydata = del_dd
dprior = gv.gvar([gv.gvar('0(1)e-04'), gv.gvar(0,L/4), gv.gvar('0(1)')])
dfit = lsq.nonlinear_fit(data=(xdata, ydata),prior=dprior,fcn=PSfcn)
print(dfit)
print("NOISY FIT")
noisydfit = lsq.nonlinear_fit(data=(xdata, ydata),prior=dprior,fcn=PSfcn, noise=True)
print("noisy reduced chi2 is %s"%(noisydfit.chi2/noisydfit.dof))

print("Fit the squared pion mass (in lattice units)")
ydata = PS_sqmass_n
mprior = gv.gvar([gv.gvar('2(2)e-03'), gv.gvar(20,10), gv.gvar('0(1)')])
mfit = lsq.nonlinear_fit(data=(xdata, ydata),prior=mprior,fcn=PSfcn)
print(mfit)

## derivatives of M^2 wrt mq eval at 1/3/5/7 ml
fitrange = np.arange(0.95,7.05,0.5)
mfitline = PSfcn(fitrange, mfit.p)
xknots = fitrange
yknots = mfitline
f = gv.cspline.CSpline(xknots, yknots)
## compute shifts
mq = [3, 5, 7]
uu = PSfcn(mq, ufit.p)
dd = PSfcn(mq, dfit.p)
delu = [ u/(i*f.D(i)) for (u,i) in zip(uu,mq) ]
deld = [ d/(i*f.D(i)) for (d,i) in zip(dd,mq) ]
print(delu)
print(deld)

chi_uu = PSfcn([1], ufit.p)[0]
chi_dd = PSfcn([1], dfit.p)[0]
chi_uu = chi_uu * overa**2
chi_dd = chi_dd * overa**2

## compute EM frac shift AT THE PHYSICAL POINT
## all quantities in physical units (GeV)
Mpi = gv.gvar('0.1349768(5)')
l3  = gv.gvar('2.81(64)')
Fpi = gv.gvar('0.1304(2)')
pi  = np.pi
D = Mpi**2 * ( 1 - (l3*Mpi**2)/16*pi**2*Fpi**2 )

del_u = chi_uu/D
del_d = chi_dd/D
print("delta_u is %s"%del_u)
print("delta_d is %s"%del_d)

print(gv.corr(del_u,del_d))
####################################################################################

amu_up =   vt['up-charge']
amu_down = vt['down-charge']

emq = gv.gvar(['3.000(1)', '5.000(1)', '7.000(1)', '27.8(1)'])
xdata = {'mq':emq}
ydata = {"up":amu_up, "down":amu_down}


prior = make_prior(xdata)
fit = lsq.nonlinear_fit(data=ydata, prior=prior, fcn=amufcn)
print("Fit results")
print(fit)

mq = np.arange(1,28,0.5)
xdata = {'mq':mq}
X = make_prior(xdata)
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
amu = amufcn(X)

## trying out splines for derivatives
yknot = amu['down']
xknot = mq
fd = gv.cspline.CSpline(xknot, yknot)
yknot = amu['up']
fu = gv.cspline.CSpline(xknot, yknot)


## 3 - 5 - 7 ml
mq = [3,5,7]

up_corr =   [ -m*d*fu.D(m) for (m,d) in zip(mq,delu) ]
down_corr = [ -m*d*fd.D(m) for (m,d) in zip(mq,deld) ]

ms = 27.8
str_corr = -ms*del_d*fd.D(ms)
print("strange correction is %s"%str_corr)
print(up_corr)
print(down_corr)

corrections = gv.BufferDict()
corrections['strange'] = str_corr
corrections['up'] = up_corr
corrections['down'] = down_corr

savefile = 'fits/em/'+a_str+'_EMcorr.p'
gv.dump(corrections, savefile, add_dependencies=True)

sys.exit(0)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

fig, ax = plt.subplots()

ax.axhline(y=0, linestyle='--')
ax.errorbar(mq, [y.mean for y in amu['down']], yerr=[y.sdev for y in amu['down']], fmt='h')
ax.errorbar(mq, [y.mean for y in amu_deriv['down']], yerr=[y.sdev for y in amu_deriv['down']], fmt='h')
plt.show()


