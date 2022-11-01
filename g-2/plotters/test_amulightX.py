import matplotlib
import sys
sys.path.append('..')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_up_diff, c_up_diff, f_up_diff
from results import vc_down_diff, c_down_diff, f_down_diff
from commonweal import delta_u, delta_d, w0overa, w0, hbarc
from gvar import log, sqrt

f_emcrc = '../fits/emcorrections/f_EMcorr.p'
c_emcrc = '../fits/emcorrections/c_EMcorr.p'
vc_emcrc = '../fits/emcorrections/vc_EMcorr.p'

f_EM_dcorr = gv.load(f_emcrc)['down']
f_EM_ucorr = gv.load(f_emcrc)['up']
c_EM_dcorr = gv.load(c_emcrc)['down']
c_EM_ucorr = gv.load(c_emcrc)['up']
vc_EM_dcorr = gv.load(vc_emcrc)['down']
vc_EM_ucorr = gv.load(vc_emcrc)['up']

up_diff = np.concatenate((vc_up_diff[:-1], c_up_diff[:-1], f_up_diff[:-1]))
down_diff = np.concatenate((vc_down_diff[:-1], c_down_diff[:-1], f_down_diff[:-1]))

print("UNCORRECTED")
print(up_diff)
print(down_diff)

vc_diff = [ (u+uEM+d+dEM) for (u,uEM,d,dEM) in zip(vc_up_diff[:-1],vc_EM_ucorr,vc_down_diff[:-1],vc_EM_dcorr) ]
c_diff = [ (u+uEM+d+dEM) for (u,uEM,d,dEM) in zip(c_up_diff[:-1],c_EM_ucorr,c_down_diff[:-1],c_EM_dcorr) ]
f_diff = [ (u+uEM+d+dEM) for (u,uEM,d,dEM) in zip(f_up_diff[:-1],f_EM_ucorr,f_down_diff[:-1],f_EM_dcorr) ]

amudiff = np.concatenate((vc_diff, c_diff, f_diff))

a = [ w0/(hbarc*w0overa[x]) for x in ['vc', 'c', 'f']]
mq = gv.gvar(['3.00(1)', '5.00(1)', '7.00(1)'])
xdata = {'a':a, 'mq':mq}

def make_prior_amudiff(x):
    prior = gv.BufferDict()
    prior['C'] = gv.gvar(['0(1)e-08', '0(1)', '0(1)'])
    prior['mq'] = x['mq'] 
    prior['a'] = x['a']
    return prior

def fcn_amudiff(p):
    out = []
    A = p['a']
    mq = p['mq']
    u0, u1, u2 = p['C']
    Y = 0.5   # Scale to make dimensionless quantity aY
    for a in A:
        for m in mq:
            #out.append( u0 * (1 + u1*(a*Y)**2 + u2*m ) )  
            out.append( u0*(1 + u1*(a*Y)**2 + u2*m) )  
    return out

prior = make_prior_amudiff(xdata)
fit = lsq.nonlinear_fit(data=amudiff, prior=prior, fcn=fcn_amudiff)
print("Fit DIFF amu results")
print(fit)

alat = [gv.gvar('0(0)')]+a
mq   = [gv.gvar('1(0)')]+[gv.gvar('3(0)')]+[gv.gvar('5(0)')]+[gv.gvar('7(0)')]
print(alat,mq)
X = make_prior_amudiff({'a':alat, 'mq':mq})
X['C'] = fit.p['C']
ext = fcn_amudiff(X)

print("a=0, mq=ml result is")
phys_diff = ext[0] 
print("up+down diff = %s"%phys_diff)

ext_a = ext[:len(mq)]      # continuum ext
ext_mq = ext[0::len(mq)]   # mq ext
print(ext_a)
print(ext_mq)
## error budget for each piece + total

#inputs = dict(up_data=up_diff, xdata=xdata)
#outputs = dict(phys_up=ext['up'][0])
#print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs, ndecimal=3, verify=True))

#inputs = dict(down_data=down_diff, xdata=xdata)
#outputs = dict(phys_down=ext['down'][0])
#print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs, ndecimal=3, verify=True))


### PLOTTING  (vs squared lattice spacing)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.gca().tick_params(right=True,top=True,direction='in')

fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, figsize=(6,3))
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
# ax1 for ml extrap
# ax2 for a^2 extrap

fitrange1 = np.arange(-0.01,1.01*a[0].mean,0.1)
X = make_prior_amudiff({'a':fitrange1, 'mq':[1]})
X['C'] = fit.p['C']
fitline1 = fcn_amudiff(X)

fitrange2 = np.arange(0.9,7.2,0.1)
X = make_prior_amudiff({'a':[0], 'mq':fitrange2})
X['C'] = fit.p['C']
fitline2 = fcn_amudiff(X)

ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext_mq], xerr=0, yerr=[y.sdev for y in ext_mq],fmt='h',mfc='none', color='blue', lw=1, marker='*')
ax2.errorbar([m.mean for m in mq], [y.mean for y in ext_a], xerr=0, yerr=[y.sdev for y in ext_a],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')

ax1.plot([p for p in fitrange1],[p.mean for p in fitline1],'--',color='r')
ax1.fill_between([p for p in fitrange1],[p.mean-p.sdev for p in fitline1],[p.mean+p.sdev for p in fitline1],alpha=0.5,lw=0,color='r')

ax2.plot([p for p in fitrange2],[p.mean for p in fitline2],'--',color='r')
ax2.fill_between([p for p in fitrange2],[p.mean-p.sdev for p in fitline2],[p.mean+p.sdev for p in fitline2],alpha=0.5,lw=0,color='r')

#ax1.legend(frameon=False, loc='center right')
#ax1.legend(frameon=False, loc='upper right')
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=4/3, color='blue', linestyle='--', lw=0.5)
ax1.set_ylabel(r'$\delta a^{(l)}_{\mu}$', rotation=0, fontsize=22, labelpad=20)
ax1.set_xlabel(r'$a^2$[GeV$^{-2}$]', fontsize=20, labelpad=10)
ax2.set_xlabel(r'$m_q/m_l$', fontsize=20, labelpad=10)
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
ax1.set_xlim(left=-0.01, right=a[0].mean*1.01)
ax2.set_xlim(left=0.9, right=7.2)
plt.tight_layout()
#plt.savefig('../figures/combX_amulightdiff.png', dpi=500)
plt.show()

plt.close()

