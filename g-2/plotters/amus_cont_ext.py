import matplotlib
import sys
sys.path.append('..')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_down_diff, c_down_diff, f_down_diff
from results import vc_amud, c_amud, f_amud
from commonweal import delta_d, hbarc, w0, hbarc, w0overa
import math as m
from gvar import log

a = [ w0/(hbarc*w0overa[x]) for x in ['vc', 'c', 'f']]
print(a)

f_emcrc = '../fits/em/f_EMcorr.p'
c_emcrc = '../fits/em/c_EMcorr.p'
vc_emcrc = '../fits/em/vc_EMcorr.p'

em_crc = gv.gvar(['-1.253(69)e-12', '-1.57(47)e-12', '-1.7(1.6)e-12'])
print(em_crc)
em_crc = gv.gvar([ gv.load(y)['strange'] for y in [vc_emcrc,c_emcrc,f_emcrc] ]) 
print(em_crc)
xdata = {'a':a}
#xdata = {'a':a, 'em':em_crc}

## directly transformed to vac pol
#vc = gv.load('../out/vc_no_vector_fit.p')
#c = gv.load('../out/c_no_vector_fit.p')
#f = gv.load('../out/f_no_vector_fit.p')
#amus_diff = [vc['down-diff']['ms'], c['down-diff']['ms'], f['down-diff']['ms']]
#amus=  [vc['down-nocharge']['ms'], c['down-nocharge']['ms'], f['down-nocharge']['ms']]

## vt fits
amus_diff_org = [vc_down_diff[-1], c_down_diff[-1], f_down_diff[-1]]
amus_diff = [(s+em) for (s,em) in zip(amus_diff_org, em_crc)]
amus      = [vc_amud[-1], c_amud[-1], f_amud[-1]]

print(amus_diff)
print(amus)

def make_prior_diff(x):
    prior = gv.BufferDict()
    prior['c'] = gv.gvar([gv.gvar(0e-12, 1e-10), gv.gvar('0(100)')])
    prior['a'] = x['a']
    #prior['em'] = x['em']
    #prior['ds'] = x['ds']
    return prior

def make_prior_tot(x):
    prior = gv.BufferDict()
    prior['c'] = gv.gvar([gv.gvar(0, 100e-10), gv.gvar('0(1)')])
    prior['a'] = x['a']
    #prior['ds'] = x['ds']
    return prior

def extrap(p):
    ans = []
    c0, c1 = p['c']
    a = p['a']
    #EM = p['em']
    #ds = p['ds']
    for p in a:
        j = c0 * ( 1 + c1*(p*0.5)**2 )
        ans.append(j)
    return ans

#amudiff
prior = make_prior_diff(xdata)
fit_diff = lsq.nonlinear_fit(data=amus_diff,prior=prior,fcn=extrap)
print(fit_diff)

contX = {'a':gv.gvar(a+['0(0)'])}
diff_ext = make_prior_diff(contX)
diff_ext['c'] = fit_diff.p['c']
diffX = extrap(diff_ext)
print("extrapolated diff is %s"%diffX)

## amus total
prior = make_prior_tot(xdata)
fit_amus = lsq.nonlinear_fit(data=amus,prior=prior,fcn=extrap)
print(fit_amus)

amus_ext = make_prior_tot(contX)
amus_ext['c'] = fit_amus.p['c']
amusX = extrap(amus_ext)
print("extrapolated amus is %s"%amusX)

## error budget
#inputs  = dict(a=a, vcdata=vc_amud, cdata=c_amud, fdata=f_amud)
#outputs = dict(amusX=amusX[0])
#print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs, ndecimal=5, verify=True))


##-----------------------------------------------------------------------#########

fitrange = np.arange(0,a[0].mean,0.001)

xx = make_prior_tot({'a':fitrange})
xx['c'] = fit_amus.p['c']
fitline_amu = extrap(xx) 

xxx = make_prior_diff({'a':fitrange})
xxx['c'] = fit_diff.p['c']
fitline_diff = extrap(xxx) 

###
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)

# diff
ax1.errorbar([i.mean**2 for i in a], [x.mean for x in amus_diff],xerr=0,yerr=[y.sdev for y in amus_diff],fmt='h',mfc='none',color='k',label='dashen corrected',lw=1, alpha=0.5)
ax1.errorbar([i.mean**2 for i in a], [x.mean for x in amus_diff_org],xerr=0,yerr=[y.sdev for y in amus_diff_org],fmt='h',mfc='none',color='r',label='uncorrected',lw=1)
ax1.plot([p**2 for p in fitrange],[p.mean for p in fitline_diff],'--',color='b')
ax1.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_diff],[p.mean+p.sdev for p in fitline_diff],alpha=0.3,lw=0,color='r')
ax1.text(x=0.35,y=1e-12, s='ext = '+str(diffX[-1]), fontsize=16)
ax1.axhline(y=0, linestyle='--')
# amus
ax2.errorbar([i.mean**2 for i in a], [x.mean for x in amus],xerr=0,yerr=[y.sdev for y in amus],fmt='h',mfc='none',color='r',label='this work',lw=1)
ax2.plot([p**2 for p in fitrange],[p.mean for p in fitline_amu],'--',color='b')
ax2.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_amu],[p.mean+p.sdev for p in fitline_amu],alpha=0.3,lw=0,color='r')
ax2.text(x=0.02,y=5.1e-9, s='ext = '+str(amusX[-1]), fontsize=16)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
ax2.legend(handles=handles,labels=labels,frameon=False,fontsize=16,loc='upper right')

ax2.set_xlabel('$a^2\ [\mathrm{GeV}^{-2}]$', fontsize=15)
ax2.set_ylabel(r'$a_{\mu}^{(s)}$', rotation=0, labelpad=20, fontsize=20)
ax1.set_ylabel(r'$\delta a_{\mu}^{(s)}$', rotation=0, labelpad=20, fontsize=20)
#plt.ylim(top=0.00)
ax2.set_xlim(left=-0.0005)

plt.tight_layout()
plt.savefig('../figures/amusX_fit.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

