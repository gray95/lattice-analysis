import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_ms_diff, c_ms_diff, f_ms_diff
from commonweal import delta_d, hbarc
import sys
import math as m
from gvar import log


a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))*1/hbarc
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))*1/hbarc
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))*1/hbarc

a = [ a_vcoarse, a_coarse, a_fine]
print(a)
ds = [ delta_d['vc'], delta_d['c'], delta_d['f'] ]
print(ds)
vc = gv.load('../out/vc_no_vector_fit.p')
c = gv.load('../out/c_no_vector_fit.p')
f= gv.load('../out/f_no_vector_fit.p')

xdata = {'a':a, 'ds':ds}

amus_diff = [vc['down-diff']['ms'], c['down-diff']['ms'], f['down-diff']['ms']]
amus_qed =  [vc['down-charge']['ms'], c['down-charge']['ms'], f['down-charge']['ms']]
amus=  [vc['down-nocharge']['ms'], c['down-nocharge']['ms'], f['down-nocharge']['ms']]
#ydiff = [vc_ms_diff, c_ms_diff, f_ms_diff]    # from vector fits

def make_prior_diff(x):
    prior = gv.BufferDict()
    prior['c'] = gv.gvar([gv.gvar(0e-12, 1e-10), gv.gvar('0(10)'), gv.gvar('0(1)')])
    prior['a'] = x['a']
    prior['ds'] = x['ds']
    return prior

def make_prior_tot(x):
    prior = gv.BufferDict()
    prior['c'] = gv.gvar([gv.gvar(0, 100e-10), gv.gvar('0(1)'), gv.gvar('0(1)')])
    prior['a'] = x['a']
    prior['ds'] = x['ds']
    return prior

def extrap(p):
    ans = []
    c0, c1, c2 = p['c']
    a = p['a']
    ds = p['ds']
    for p,q in zip(a,ds):
        j = c0 * ( 1 + c1*(p*0.5/m.pi)**2 + c2*(q/(1-q)) )
        ans.append(j)
    return ans
#amudiff
prior = make_prior_diff(xdata)
fit = lsq.nonlinear_fit(data=amus_diff,prior=prior,fcn=extrap)
print(fit)
fit.p['a'] = 0
fit.p['ds'] = delta_d['c']
diff_ext = extrap(fit.p)
print(diff_ext)
## amu
prior = make_prior_tot(xdata)
fit = lsq.nonlinear_fit(data=amus_qed,prior=prior,fcn=extrap)
print(fit)
prior = make_prior_tot(xdata)
fit = lsq.nonlinear_fit(data=amus,prior=prior,fcn=extrap)
print(fit)

## error budget
inputs  = dict(a=a, ds=ds, diff=amus_diff)
outputs = dict(ext=diff_ext[0])
print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs))

##-----------------------------------------------------------------------#########

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,a_vcoarse.mean,0.001)

xx = make_prior_tot({'a':fitrange, 'ds':delta_d['c']*np.ones(len(fitrange))})
xx['c'] = fit.p['c']
fitline = extrap(xx) 

Y = make_prior_tot({'a':[0], 'ds':[delta_d['c'].mean]})
Y['c'] = fit.p['c']
#Y['C'] = fit.p['C']
plym_rt = extrap(Y)
print("continuum extrapolation: ", plym_rt)
plt.errorbar([i.mean**2 for i in a], [x.mean for x in amus],xerr=0,yerr=[y.sdev for y in amus],fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='b')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.3,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{GeV}^{-2}]$', fontsize=15)
plt.ylabel(r'$a_{\mu}^{(s)}$', rotation=0, labelpad=20, fontsize=20)
#plt.title('Preliminary continuum extrapolation')
#plt.ylim(top=0.00)
plt.xlim(left=-0.0005)

plt.tight_layout()
plt.savefig('../figures/amusX.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

