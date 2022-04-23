import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_up_diff, c_up_diff, f_up_diff
from results import vc_down_diff, c_down_diff, f_down_diff
from commonweal import delta_u, delta_d, w0overa, w0, hbarc
from gvar import log
import sys

#vc = gv.load('../out/vc_no_vector_fit.p')
#c = gv.load('../out/c_no_vector_fit.p')
#f= gv.load('../out/f_no_vector_fit.p')

#keys = ['3ml','5ml','7ml']
#vc_up_diff = [vc['up-diff'][k] for k in keys]
#vc_down_diff = [vc['down-diff'][k] for k in keys]
#c_up_diff = [c['up-diff'][k] for k in keys]
#c_down_diff = [c['down-diff'][k] for k in keys]
#f_up_diff = [f['up-diff'][k] for k in keys]
#f_down_diff = [f['down-diff'][k] for k in keys]
  
#up_diff = np.concatenate((vc_up_diff, c_up_diff, f_up_diff))
#down_diff = np.concatenate((vc_down_diff, c_down_diff, f_down_diff))
## correlated vector fits

up_diff = np.concatenate((vc_up_diff, c_up_diff, f_up_diff))
down_diff = np.concatenate((vc_down_diff, c_down_diff, f_down_diff))
print(up_diff)
print(down_diff)

a = [ w0/(hbarc*w0overa[x]) for x in ['vc', 'c', 'f']]
du = [delta_u['vc'], delta_u['c'], delta_u['f']] 
dd = [delta_d['vc'], delta_d['c'], delta_d['f']] 

emq = gv.gvar(['3.000(1)', '5.000(1)', '7.000(1)'])
#emq = [3,5,7]
xdata = {'a':a, 'mq':emq, 'dd':dd, 'du':du}
ydata = {"up":up_diff, "down":down_diff}


def make_prior(x):
    prior = gv.BufferDict()
    prior['uC'] = gv.gvar(['0(1)e-9', '0(5)', '0(5)'])
    prior['dC'] = gv.gvar(['0(1)e-9', '0(5)', '0(5)'])
    prior['mq'] = x['mq'] 
    prior['a'] = x['a']
    prior['dd']= x['dd']
    prior['du']= x['du']
    return prior

def fcn(p):
    out = {'up':[], 'down':[]}
    du = p['du']
    dd = p['dd']
    A = p['a']
    M = p['mq']
    for (a,d,u) in zip(A,dd,du):
        for m in M:
            out['up'].append( p['uC'][0] * (1 + p['uC'][1]*m + p['uC'][2]*a**2 ) )  
            out['down'].append( p['dC'][0] * (1 + p['dC'][1]*m + p['dC'][2]*a**2) ) 
    return out

prior = make_prior(xdata)
fit = lsq.nonlinear_fit(data=ydata, prior=prior, fcn=fcn)
print("Fit results")
print(fit)

alat = [gv.gvar('0(0)')]+a
mq   = [1, 3, 5, 7]
print(alat)
X = make_prior({'a':alat, 'mq':mq,
                'dd':len(alat)*[delta_d['c']], 
                'du':len(alat)*[delta_u['c']]})

X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
ext = fcn(X)
print("a=0, mq=ml result is")
phys_diff = ext['down'][0] + ext['up'][0]
print(phys_diff)
ext_a = {'up':ext['up'][::4], 'down':ext['down'][::4]}      # continuum ext
ext_mq = {'up':ext['up'][:4], 'down':ext['down'][:4]}    # mq ext
## ETM number
## renorm QED+QCD quark mass == renorm QCD quark mass
##    => xt to mq=ml and mq=ml/Z_m^QED


### PLOTTING  (vs squared lattice spacing)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

fig, (ax1, ax2) = plt.subplots(1,2,sharey=True)
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
# ax1 for ml extrap
# ax2 for a^2 extrap

fitrange1 = np.arange(-0.01,1.01*a[0].mean,0.01)
X = make_prior({'a':fitrange1, 'mq':[1], 'dd':delta_d['c']*np.ones(len(fitrange1)), 'du':delta_u['c']*np.ones(len(fitrange1))})
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
fitline1 = fcn(X)

fitrange2 = np.arange(0.9,7.2,0.01)
X = make_prior({'a':[0], 'mq':fitrange2, 'dd':delta_d['c']*np.ones(len(fitrange2)), 'du':delta_u['c']*np.ones(len(fitrange2))})
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
fitline2 = fcn(X)

ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext_mq['down']], xerr=0, yerr=[y.sdev for y in ext_mq['down']],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')
ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext_mq['up']], xerr=0, yerr=[y.sdev for y in ext_mq['up']],fmt='h',mfc='none', label='up like', color='red', lw=1, marker='*')


ax2.errorbar(mq, [y.mean for y in ext_a['down']], xerr=0, yerr=[y.sdev for y in ext_a['down']],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')
ax2.errorbar(mq, [y.mean for y in ext_a['up']], xerr=0, yerr=[y.sdev for y in ext_a['up']],fmt='h',mfc='none', label='up like', color='red', lw=1, marker='*')

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
ax1.set_ylabel(r'$\delta a^{(l)}_{\mu}$', rotation=0, fontsize=25, labelpad=20)
ax1.set_xlabel(r'$a^2$', fontsize=25, labelpad=10)
ax2.set_xlabel(r'$m_q/m_l$', fontsize=25, labelpad=10)
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
ax1.set_xlim(left=-0.01, right=a[0].mean*1.01)
ax2.set_xlim(left=0.9, right=7.2)
plt.tight_layout()
plt.savefig('../figures/TESTcombX_amulightdiff.png', dpi=500)
plt.show()

plt.close()

