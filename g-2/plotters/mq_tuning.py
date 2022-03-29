
import matplotlib.pyplot as plt
import numpy as np
import gvar as gv
from scipy.optimize import curve_fit
import math
import lsqfit as lsq
import sys
from commonweal import w0, w0overa, hbarc

#M_etas     = gv.gvar('0.6885(22)')        # GeV
#M_etas_qed = gv.gvar('0.68989(49)')     # BMW 
#aM_etas = M_etas * (1/hbarc) * (w0/w0overa['f'])
#print(aM_etas**2)
#am = [0.0362, 0.0364, 0.0366, 0.0368]
am = [3, 5, 7]#, 27.5]
overa = (w0overa['vc']/w0)*hbarc

obs = gv.load('../pseudoscalar/ps_vcoarse.p')
#obs = gv.load('../pseudoscalar/ps_coarse.p')
#ratio = [ a**2/b**2 for (a,b) in zip(obs['Eqed'],obs['E']) ]
#print(ratio)
PS_sqmass_n = [ y**2 for y in obs['E:n'] ]
PS_sqmass_d = [ y**2 for y in obs['E:d'] ]
PS_sqmass_u = [ y**2 for y in obs['E:u'] ]
del PS_sqmass_n[-1]
del PS_sqmass_u[-1]
del PS_sqmass_d[-1]

del_uu = [ (A**2-C**2)*overa**2 for (A,C) in zip(obs['E:u'],obs['E:n']) ]
del_dd = [ (B**2-C**2)*overa**2 for (B,C) in zip(obs['E:d'],obs['E:n']) ]
del del_uu[-1]
del del_dd[-1]

xdata = am #np.array([ x for x in am ])
ydata = del_uu
y_data = [ y.mean for y in del_uu ]
y_err  = [ y.sdev for y in del_uu ]

#tuned_eta = aM['e/3'][1]['E'] * hbarc * (w0overa['f']/w0) 
#print("Compare %s with %s" % (M_etas*1000, tuned_eta*1000))
#ms_qed = 0.0364 * (ratio[1])**(-2)
#print(ms_qed)
#print((ms_qed-0.0364)/0.0364)

##########################################################################
#Fitting

def extrap(x,p):
    res = []
    for pt in x:
        res.append( p[0]*(1 + p[1]*pt + p[2]*pt**2 ) )
    return res

def fitargs(z):
    p0  = z[0]
    dp0 = z[1]
    p1  = z[2]
    dp1 = z[3]
    prior = gv.gvar([gv.gvar(p0, dp0), gv.gvar(p1, dp1)])
    return dict(prior=prior, fcn=extrap, data=(xdata, ydata))

##########################################################################

#ydata = PS_sqmass
#z0 = [1.0, 0.3, 0, 0.1] 
#fit, z = lsq.empbayes_fit(z0, fitargs)
#print(fit)
#tuned_ms = (1-( (fit.p[0]+fit.p[1]*0.0364) -1))*0.0364
#print(tuned_ms)
#print(fit.p[0]-1)
#print("prior width that maxmises logGBF: ", z)
## renormalise with QED
#Zm_qed = gv.gvar('1.000724(11)')
#print( "quark mass WITH QED is %s MeV" % (overa*tuned_ms*1000*Zm*Zm_qed) )

# with quark mass squared term - NO emp bayes
print("Fitting UP squared diff")
print("------------------------------------------------------")
prior = {}
prior = gv.gvar([gv.gvar('0.0001(2)'), gv.gvar('10(50)'), gv.gvar('0(1)')])
ufit = lsq.nonlinear_fit(data=(xdata, ydata),prior=prior,fcn=extrap, debug=True)
#print("WITH SQUARED QUARK MASS TERM\n")
print(ufit)

print("Fitting DOWN squared diff")
print("------------------------------------------------------")
ydata = del_dd
prior = gv.gvar([gv.gvar('0.0001(2)'), gv.gvar('10(50)'), gv.gvar('0(1)')])
dfit = lsq.nonlinear_fit(data=(xdata, ydata),prior=prior,fcn=extrap, debug=True)
#print("WITH SQUARED QUARK MASS TERM\n")
print(dfit)

xx = np.linspace(xdata[-1], 1, 100)
fitrange = np.linspace(xdata[-1], 1, 100)

ufitline = extrap(fitrange,ufit.p)
dfitline = extrap(fitrange,dfit.p)

chi_uu = extrap([1], ufit.p)[0]
chi_dd = extrap([1], dfit.p)[0]
print("extrapolated value at mq/ml = %f is %s GeV^2\n(not renormalised)"%(1,chi_uu))
print("extrapolated value at mq/ml = %f is %s GeV^2\n(not renormalised)"%(1,chi_dd))

## compute EM frac shift
## all quantities in physcial units (GeV)
Mpi = gv.gvar('0.1349768(5)')   
l3  = gv.gvar('2.81(64)')
Fpi = gv.gvar('0.1304(2)')
pi  = math.pi
D = Mpi**2 * ( 1 - (l3*Mpi**2)/16*pi**2*Fpi**2 )  

del_u = chi_uu/D
del_d = chi_dd/D
print(del_u)
print(del_d)
#######PLOTTING

plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['New Century Schoolbook'],'size':12})
plt.rc('axes', linewidth=0.5)

fig, ax1 = plt.subplots()
#fig, (ax1, ax2) = plt.subplots( 2, 1, sharex=True )
#ax1.plot(xx, f_model, "r--", label='linear fit')
#ax1.set_xlim(left=0.99*xdata[0], right=1.01*xdata[-1])
ax1.set_xlabel('$m_q/m_l$', fontsize=20)
#ax1.set_ylabel(r'$R^0_{\mbox{\scriptsize QED}}[M^2_{\eta_s}]$', rotation=0, labelpad=40, fontsize=20)
ax1.set_ylabel('$\Delta M_{xx}^2$\n[GeV]$^2$', rotation=0, labelpad=30, fontsize=20)
ax1.errorbar( [x for x in xdata], [y.mean for y in del_uu], yerr=[y.sdev for y in del_uu], color='red', fmt='h', label=r'$\Delta Q=2e/3$')
ax1.errorbar( [x for x in xdata], [y.mean for y in del_dd], yerr=[y.sdev for y in del_dd], color='blue', fmt='h', label=r'$\Delta Q=e/3$')


ax1.fill_between([x for x in fitrange], [y.mean-y.sdev for y in dfitline], [y.mean+y.sdev for y in dfitline], alpha=0.1,color='blue')
ax1.fill_between([x for x in fitrange], [y.mean-y.sdev for y in ufitline], [y.mean+y.sdev for y in ufitline], alpha=0.1,color='red')
ax1.legend(frameon=False, loc='upper left')

plt.title('0.15fm - 1844 cfgs')
plt.tight_layout()
#plt.savefig('../figures/mq_emtuning_fine.png', dpi=500)
plt.show()

