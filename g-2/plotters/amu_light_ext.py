import matplotlib
import sys
sys.path.append('..')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_up_diff, c_up_diff, f_up_diff
from results import vc_down_diff, c_down_diff, f_down_diff
from commonweal import w0overa, w0, hbarc
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

# correct for dashen scheme
vc_up_diff = [ (y+corr) for (y,corr) in zip(vc_up_diff[:-1],vc_EM_ucorr) ]
vc_down_diff = [ (y+corr) for (y,corr) in zip(vc_down_diff[:-1],vc_EM_dcorr) ]

c_up_diff = [ (y+corr) for (y,corr) in zip(c_up_diff[:-1],c_EM_ucorr) ]
c_down_diff = [ (y+corr) for (y,corr) in zip(c_down_diff[:-1],c_EM_dcorr) ]

f_up_diff = [ (y+corr) for (y,corr) in zip(f_up_diff[:-1],f_EM_ucorr) ]
f_down_diff = [ (y+corr) for (y,corr) in zip(f_down_diff[:-1],f_EM_dcorr) ]

up_diff = np.concatenate((vc_up_diff, c_up_diff, f_up_diff))
down_diff = np.concatenate((vc_down_diff, c_down_diff, f_down_diff))

print(up_diff)
print(down_diff)

ENSEMBLES = ['vc', 'c', 'f']
a = [ w0/(hbarc*w0overa[x]) for x in ENSEMBLES]

emq = gv.gvar(['3.0(1)', '5.0(1)', '7.0(1)'])#, '27.8(1)'])
xdata = {'a':a, 'mq':emq}
amudiff_data = {"up":up_diff, "down":down_diff}

def make_prior_amudiff(x):
    prior = gv.BufferDict()
    prior['uC'] = gv.gvar(['0(1)e-9', '0(1)', '0(1)'])
    prior['dC'] = gv.gvar(['0(1)e-10', '0(1)', '0(1)'])
    prior['mq'] = x['mq'] 
    prior['a'] = x['a']
    return prior

def fcn_amudiff(p):
    out = {'up':[], 'down':[]}
    A = p['a']
    M = p['mq']
    u0, u1, u2 = p['uC']
    d0, d1, d2 = p['dC']
    Y = 0.5  # Scale to make dimensionless quantity aY
    for a in A:
        for m in M:
            out['up'].append(   u0*(1 + u1*(a*Y)**2 + u2*m  ))   
            out['down'].append( d0*(1 + d1*(a*Y)**2 + d2*m  )) 
    return out


prior = make_prior_amudiff(xdata)
fit = lsq.nonlinear_fit(data=amudiff_data, prior=prior, fcn=fcn_amudiff)
print("Fit DIFF amu results")
print(fit)

alat = [gv.gvar('0(0)')]+a
mq   = gv.gvar(['1(0)', '3(0)', '5(0)', '7(0)'])#, '27.8(0)'])
X = make_prior_amudiff({'a':alat, 'mq':mq})

X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
ext = fcn_amudiff(X)

print("a=0, mq=ml result is")
print("down diff = %s"%ext['down'][0])
print("up diff = %s"%ext['up'][0])
phys_diff = ext['down'][0] + ext['up'][0]
print("up+down diff = %s"%phys_diff)

ext_a = {'up':ext['up'][:len(mq)], 'down':ext['down'][:len(mq)]}      # continuum ext
ext_mq = {'up':ext['up'][0::len(mq)], 'down':ext['down'][0::len(mq)]}    # mq ext
print(ext_a)
print(ext_mq)
sys.exit(0)
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
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
fitline1 = fcn_amudiff(X)

fitrange2 = np.arange(0.9,7.2,0.1)
X = make_prior_amudiff({'a':[0], 'mq':fitrange2})
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
fitline2 = fcn_amudiff(X)

ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext_mq['down']], xerr=0, yerr=[y.sdev for y in ext_mq['down']],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')
ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext_mq['up']], xerr=0, yerr=[y.sdev for y in ext_mq['up']],fmt='h',mfc='none', label='up like', color='red', lw=1, marker='*')

ax2.errorbar([m.mean for m in mq], [y.mean for y in ext_a['down']], xerr=0, yerr=[y.sdev for y in ext_a['down']],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')
ax2.errorbar([m.mean for m in mq], [y.mean for y in ext_a['up']], xerr=0, yerr=[y.sdev for y in ext_a['up']],fmt='h',mfc='none', label='up like', color='red', lw=1, marker='*')

ax1.plot([p for p in fitrange1],[p.mean for p in fitline1['up']],'--',color='r')
ax1.fill_between([p for p in fitrange1],[p.mean-p.sdev for p in fitline1['up']],[p.mean+p.sdev for p in fitline1['up']],alpha=0.5,lw=0,color='r')
ax1.plot([p for p in fitrange1],[p.mean for p in fitline1['down']],'--',color='blue')
ax1.fill_between([p for p in fitrange1],[p.mean-p.sdev for p in fitline1['down']],[p.mean+p.sdev for p in fitline1['down']],alpha=0.5,lw=0,color='blue')

ax2.plot([p for p in fitrange2],[p.mean for p in fitline2['up']],'--',color='r')
ax2.fill_between([p for p in fitrange2],[p.mean-p.sdev for p in fitline2['up']],[p.mean+p.sdev for p in fitline2['up']],alpha=0.5,lw=0,color='r')
ax2.plot([p for p in fitrange2],[p.mean for p in fitline2['down']],'--',color='blue')
ax2.fill_between([p for p in fitrange2],[p.mean-p.sdev for p in fitline2['down']],[p.mean+p.sdev for p in fitline2['down']],alpha=0.5,lw=0,color='blue')

#ax1.legend(frameon=False, loc='center right')
#ax1.legend(frameon=False, loc='upper right')
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=4/3, color='blue', linestyle='--', lw=0.5)
ax1.set_ylabel(r'$\delta a^{(l)}_{\mu}$', rotation=0, fontsize=20, labelpad=20)
ax1.set_xlabel(r'$a^2$[GeV$^{-2}$]', fontsize=20, labelpad=10)
ax2.set_xlabel(r'$m_q/m_l$', fontsize=20, labelpad=10)
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
ax1.set_xlim(left=-0.01, right=a[0].mean*1.01)
ax2.set_xlim(left=0.9, right=7.2)
plt.tight_layout()
plt.savefig('../figures/combX_amulightdiff.png', dpi=500)
plt.show()

plt.close()

