import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from commonweal import bmw_diff
from results import vc_diff_phys_ext, c_diff_phys_ext
from results import vc_amus_diff, c_amus_diff, f_amus_diff
import sys
import math as m
from gvar import log

bmw = bmw_diff*10**-10

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))

a = [ a_vcoarse, a_coarse, a_fine]
print(a)
xdata = np.array([x.mean for x in a])

ydata = [vc_diff_phys_ext, c_diff_phys_ext, 0.8*c_diff_phys_ext]

# correlation
#vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1]])
####--------------------------------------------------------------------#########
#def extrap(x,p):
#    res = []
#    for i,point in enumerate(x):
#        res.append(p['a'][0] * (1 + p['a'][1]*(point*0.5)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
#    return res

def extrap(x,p):
    return ( p[0] + p[1]*(x)**2 ) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))

def fitargs(z):
    dp = z
    prior = gv.gvar([gv.gvar(-1e-11, dp), gv.gvar(0, dp)])
    return dict(prior=prior, fcn=extrap, data=(xdata, ydata))

##-----------------------------------------------------------------------#########

## EMPIRICAL BAYES
#z0 = { 'dp0':10, 'dp1':10 }
z0 = 1e-6
fit, z = lsq.empbayes_fit(z0, fitargs)
print(fit.format(True))
print("prior on a[0] that maxmises logGBF: ", z)
## No emp bayes
#prior = {}
#prior['a'] = [gv.gvar('1(0.5)'),gv.gvar(0,0.1)]
#fit = lsq.nonlinear_fit(data=(xdata, ydata),prior=prior,fcn=extrap)
#print(fit)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,a_vcoarse.mean,0.00001)
fitline = [extrap(x, fit.p) for x in fitrange]

plym_rt = extrap(0, fit.p)
#print("BMW " , bmw)
print("continuum extrapolation: ", plym_rt)
#print("our error / bmw error: ", plym_rt[0].sdev/bmw.sdev)

plt.errorbar([i.mean**2 for i in a], [x.mean for x in ydata],xerr=0,yerr=[y.sdev for y in ydata],fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='b')
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

