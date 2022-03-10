import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from commonweal import bmw_diff, hbarc
#from results import vc_diff_phys_ext, c_diff_phys_ext
from results import vc_diff, c_diff
import sys
import math as m
from gvar import log

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))*1/hbarc
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))*1/hbarc
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))*1/hbarc

a = [a_vcoarse, a_coarse, a_fine]
mq = gv.gvar(['3(0.001)', '5(0.001)', '7(0.001)'])
#msqed = [gv.gvar('-0.0014(1)'), gv.gvar('-0.001314(83)'), gv.gvar('-0.001249(35)')]
xdata = {'a':a, 'mq':mq}
f_diff = [gv.gvar('-4.0(9.4)e-11'), gv.gvar('-6.1(8.2)e-11'), gv.gvar('-5.2(4.1)e-11')]
diff = np.concatenate((vc_diff, c_diff, f_diff))
print(diff)
ydata = {'y1':diff}

def make_prior(x):
    prior = gv.BufferDict()
    prior['c'] = gv.gvar([gv.gvar(0, 1e-9), gv.gvar(0, 1), gv.gvar(0, 1)])
    prior['a'] = x['a']
    prior['mq'] = x['mq']
    return prior

def extrap(p):
    ans = {'y1':[]}
    c0, c1, c2 = p['c']
    a = p['a']
    mq = p['mq']
    for i in a:
        j = c0 * ( 1 + c1*(i*0.5)**2 + c2*(mq) ) 
        ans['y1'].append(j)
    return ans 

def fitargs(z):
    dp1 = 1e-8
    dp2 = z[0]
    dp3 = z[1]
    prior = make_prior(xdata)
    prior['c'] = gv.gvar([gv.gvar(0, dp1), gv.gvar(0, dp2), gv.gvar(0, dp3)])
    return dict(prior=prior, fcn=extrap, data=ydata)

##-----------------------------------------------------------------------#########

## EMPIRICAL BAYES
#z0 = [10, 10]
#z0 = { 'param':1e3, "data":11.1 }
#fit, z = lsq.empbayes_fit(z0, fitargs)
#print(fit.format(True))
#print("prior on a[0] that maxmises logGBF: ", z)

## No emp bayes
prior = make_prior(xdata)
fit = lsq.nonlinear_fit(data=ydata,prior=prior,fcn=extrap)


print(fit)
print('\n================ Add noise to prior')
noisyfit = lsq.nonlinear_fit( data=ydata, prior=prior, fcn=extrap, noise=True)
print(noisyfit)
####################
fit.p['a'] = [0]
fit.p['mq'] = 1
physical_result= fit.fcn(fit.p)
noisyfit.p['a'] = [0]
noisyfit.p['mq'] = 1
noisy_phys_res= noisyfit.fcn(noisyfit.p)
print("physical point result is %s"%physical_result)
print("physical point result is %s"%noisy_phys_res)
sys.exit(0)

#### PLOT SOME TASTY FIGS
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,a_vcoarse.mean,0.0001)
fitline = [ fit.p['c'][0]*(1+fit.p['c'][1]*x**2+fit.p['c'][2]*msqed[2]) for x in fitrange ]

plt.errorbar([i.mean**2 for i in a], [y.mean for y in ydata],xerr=0,yerr=[y.sdev for y in ydata],fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.3,lw=0,color='r')

#plt.errorbar(0,bmw.mean,xerr=0,yerr=bmw.sdev,fmt='s',color='b',label='BMW',lw=1)
#plt.errorbar(0,giu_diff.mean,xerr=0,yerr=giu_diff.sdev,fmt='s',color='g',label='Giusti et al.',lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.ylabel(r'$\delta a_{\mu}^{\small\mathrm{qed,conn}}$', rotation=0, labelpad=20, fontsize=20)
#plt.title('Preliminary continuum extrapolation')
#plt.ylim(top=0.00)
plt.xlim(left=-0.0005)

plt.tight_layout()
#plt.savefig('../figures/MADEUPdiff_continuum_ext.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

