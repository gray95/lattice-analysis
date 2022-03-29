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

up_diff = np.concatenate((vc_up_diff, c_up_diff, f_up_diff))
down_diff = np.concatenate((vc_down_diff, c_down_diff, f_down_diff))
a = [ w0/(hbarc*w0overa[x]) for x in ['vc', 'c', 'f']]
du = [delta_u['vc'], delta_u['c'], delta_u['f']] 
dd = [delta_d['vc'], delta_d['c'], delta_d['f']] 

emq = gv.gvar(['3.000(1)', '5.000(1)', '7.000(1)'])
xdata = {'a':a, 'mq':emq, 'dd':dd, 'du':du}
ydata = {"up": up_diff, "down":down_diff}


def make_prior(x):
    prior = gv.BufferDict()
    prior['uC'] = gv.gvar(['0(1)e-9', '1(2)', '1(2)'])
    prior['dC'] = gv.gvar(['0(1)e-9', '1(2)', '1(2)'])
    prior['mq']= x['mq']
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
            out['up'].append( p['uC'][0] * (1 + p['uC'][1]*m*u + p['uC'][2]*a**2 ) )  
            out['down'].append( p['dC'][0] * (1 + p['dC'][1]*m*d + p['dC'][2]*a**2) ) 
    return out

prior = make_prior(xdata)
fit = lsq.nonlinear_fit(data=ydata, prior=prior, fcn=fcn)
print("Fit results")
print(fit)

alat = [gv.gvar('0(0)')]+a
print(alat)
X = make_prior({'a':alat, 'mq':[1],
                'dd':len(alat)*[delta_d['c']], 
                'du':len(alat)*[delta_u['c']]})

X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
ext = fcn(X)

### PLOTTING  (vs squared lattice spacing)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

#fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
fig, ax1 = plt.subplots()

fitrange = np.arange(-0.01,a[0].mean,0.001)
X = make_prior({'a':fitrange, 'mq':[1], 'dd':delta_d['c']*np.ones(len(fitrange)), 'du':delta_u['c']*np.ones(len(fitrange))})
X['uC'] = fit.p['uC']
X['dC'] = fit.p['dC']
fitline = fcn(X)

#ax1.errorbar([x for x in xdata] ,[y.mean for y in ydata["up"]],xerr=0,yerr=[y.sdev for y in ydata["up"]] ,fmt='h',color='r',lw=1, label="$Q=2e/3$")
#ax1.errorbar([x.mean for x in xdata['a']] ,[y.mean for y in ydata["down"]],xerr=0,yerr=[y.sdev for y in ydata["down"]] ,fmt='h',color='blue',lw=1, label="$Q=e/3$")

# extrapolated 
ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext['down']], xerr=[x.sdev**2 for x in alat], yerr=[y.sdev for y in ext['down']],fmt='h',mfc='none', label='down like', color='blue', lw=1, marker='*')
ax1.errorbar([x.mean**2 for x in alat], [y.mean for y in ext['up']], xerr=[x.sdev**2 for x in alat], yerr=[y.sdev for y in ext['up']],fmt='h',mfc='none', label='up like', color='red', lw=1, marker='*')

ax1.plot([p for p in fitrange],[p.mean for p in fitline['up']],'--',color='r')
ax1.fill_between([p for p in fitrange],[p.mean-p.sdev for p in fitline['up']],[p.mean+p.sdev for p in fitline['up']],alpha=0.5,lw=0,color='r')
ax1.plot([p for p in fitrange],[p.mean for p in fitline['down']],'--',color='blue')
ax1.fill_between([p for p in fitrange],[p.mean-p.sdev for p in fitline['down']],[p.mean+p.sdev for p in fitline['down']],alpha=0.5,lw=0,color='blue')

ax1.legend(frameon=False, loc='center right')
#ax1.legend(frameon=False, loc='upper right')
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=4/3, color='blue', linestyle='--', lw=0.5)

ax1.set_ylabel(r'$\delta a^{(l)}_{\mu}$', rotation=0, fontsize=25, labelpad=20)
ax1.set_xlabel(r'$a^2$', fontsize=25, labelpad=10)
#plt.ylabel(r'$\frac{a_{\mu}^{\mathrm{qcd+qed}}}{a_{\mu}^{\mathrm{qcd}}}$', rotation=0, labelpad=25, fontsize=20)
#plt.ylabel(r'$a^2M^2_{\eta}$', rotation=0, labelpad=25, fontsize=20)
ax1.set_title('Continuum Extrapolation with $m_q=m_l$')
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
ax1.set_xlim(left=-0.01)

plt.tight_layout()
plt.savefig('../figures/comb_ext_amulight.png', dpi=500)
plt.show()

plt.close()

