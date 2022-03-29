import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from results import vc_up_noqed, vc_down_noqed, vc_amu_up, vc_amu_down
from results import vc_up_phys_diff, vc_down_phys_diff, vc_up_diff, vc_down_diff
from results import vc_phys15
from results import vc_amu 
from commonweal import w0, w0overa, hbarc
from gvar import log
import sys

M_etas     = gv.gvar('0.6885(22)')        # GeV
M_etas_qed = gv.gvar('0.68989(49)')     # BMW 
overa = (w0overa['f']/w0)*hbarc

aM_ss = M_etas * 1/overa

print(aM_ss**2)

ms_tuning = gv.load('../pseudoscalar/mistune_ps_fine.p')
amus_diff = [gv.gvar('-4.42(34)e-12'),gv.gvar('-4.53(36)e-12'),gv.gvar('-4.40(34)e-12'),gv.gvar('-4.39(34)e-12')]   # correlated fit
#xdata = np.array( [3, 5, 7] )
xdata = np.array( [0.0362, 0.0364, 0.0366, 0.0368] )
#ydata = {"up": vc_amu_up, "down":vc_amu_down}
ydata = {   "Q=0":[y**2 for y in ms_tuning['E:n']], 
        "Q=0.101":[y**2 for y in ms_tuning['E:d']], 
        "Q=0.202":[y**2 for y in ms_tuning['E:u']]  }

def fcn(x,p):
    out = {'Q=0':[], 'Q=0.101':[], 'Q=0.202':[]}
    for pt in x:
        out['Q=0'].append( p['Q=0'][0] + p['Q=0'][1]*pt )  
        out['Q=0.101'].append(p['Q=0.101'][0] + p['Q=0.101'][1]*pt )  
        out['Q=0.202'].append(p['Q=0.202'][0] + p['Q=0.202'][1]*pt )  
    return out

def fitargs(z):
    dp1 = z[0]
    dp2 = z[1]
    prior = {}
    prior['up'] = [gv.gvar(0, 1e-9), gv.gvar(0, 1)]
    prior['down'] = [gv.gvar(0, 1e-9), gv.gvar(0, 1)]
    return dict(prior=prior, fcn=fcn, data=(xdata, ydata))

## EMPIRICAL BAYES
#z0 = { 'dp0' : 0.1, 'dp1' : '0.0' }
#z0 = [10, 1]
#fit, z = lsq.empbayes_fit(z0, fitargs)
#print(fit.format(True))
#print("prior width that maxmises logGBF: ", z)
#print("prior width on a[1] that maxmises logGBF: ", z['dp1'])

prior = {}
prior['Q=0'] = [gv.gvar('0(1)'),gv.gvar('3(3)')]
prior['Q=0.101'] = [gv.gvar('0(1)'),gv.gvar('3(3)')] 
prior['Q=0.202'] = [gv.gvar('0(1)'),gv.gvar('3(3)')] 
fit = lsq.nonlinear_fit(data=(xdata, ydata), prior=prior, fcn=fcn)
print("Fit results")
print(fit)

print(fit.p)
tuned_ms = (aM_ss**2 - fit.p['Q=0.101'][0])*1/fit.p['Q=0.101'][1]
print("tuned quark mass is %s"%tuned_ms)

##########################################################################
############# PLOTTING ###################################################
##########################################################################
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
#plt.gca().tick_params(right=True,top=True,direction='in')

#fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
fig, ax1 = plt.subplots()

fitrange = np.arange(xdata[0],xdata[-1],1e-6)
fitline = fcn(fitrange ,fit.p)

ax1.errorbar( xdata ,[y.mean for y in ydata["Q=0"]],xerr=0,yerr=[y.sdev for y in ydata["Q=0"]] ,fmt='h',color='black',lw=1, label="$Q=0$")
ax1.plot([p for p in fitrange],[p.mean for p in fitline['Q=0']], '--',color='black')
ax1.errorbar( xdata ,[y.mean for y in ydata["Q=0.101"]],xerr=0,yerr=[y.sdev for y in ydata["Q=0.101"]] ,fmt='h',color='blue',lw=1, label="$Q=-e/3$")
ax1.plot([p for p in fitrange],[p.mean for p in fitline['Q=0.101']], '--',color='blue')
ax1.errorbar( xdata ,[y.mean for y in ydata["Q=0.202"]],xerr=0,yerr=[y.sdev for y in ydata["Q=0.202"]] ,fmt='h',color='red',lw=1, label="$Q=0$")
ax1.plot([p for p in fitrange],[p.mean for p in fitline['Q=0.202']], '--',color='red')

# extrapolated 
#ax1.errorbar(2/3+0.01, ext_up.mean,xerr=0,yerr=ext_up.sdev,fmt='h',mfc='none', label='$m_u$ linear extrapolation', color='red', lw=1, marker='*')

ax1.fill_between( fitrange, [p.mean-p.sdev for p in fitline['Q=0.202']],[p.mean+p.sdev for p in fitline['Q=0.202']],alpha=0.5,lw=0,color='red')
#ax1.plot([p for p in fitrange],[p['down'].mean for p in fitline],'--',color='blue')

## phys values
#ax1.errorbar(2/3-0.02, vc_up_phys_diff.mean,xerr=0,yerr=vc_up_phys_diff.sdev,fmt='h',label='$m_u, m_d$ sim $t^* = 1.5$fm', color='black', lw=1, marker='s', markersize=4)
#ax1.errorbar(4/3-0.02, vc_down_phys_diff.mean,xerr=0,yerr=vc_down_phys_diff.sdev,fmt='h',color='black', lw=1, marker='s', markersize=4)

ax1.legend(frameon=False, loc='lower right')
#ax1.legend(frameon=False, loc='upper right')
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=4/3, color='blue', linestyle='--', lw=0.5)

#ax1.set_xticks([2/3,4/3,3,5,7])
#ax1.set_xticklabels(['$m_u$','$m_d$','3','5','7'])
#ax1.set_ylabel(r'$\delta a^{\tiny\mbox{HVP}}_{\mu}$', rotation=0, fontsize=15, labelpad=10)
ax1.set_xlabel(r'$\frac{m_q}{m_l}$', fontsize=25, labelpad=30)
#ax2.set_xlabel(r'$am_q$', fontsize=25, labelpad=20)
#plt.ylabel(r'$\frac{a_{\mu}^{\mathrm{qcd+qed}}}{a_{\mu}^{\mathrm{qcd}}}$', rotation=0, labelpad=25, fontsize=20)
#plt.ylabel(r'$a^2M^2_{\eta}$', rotation=0, labelpad=25, fontsize=20)
ax1.set_title('0.09fm - 50 cfgs - correlated fit')
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
#ax1.set_ylim(bottom=-3e-10, top=0.8e-10)

plt.tight_layout()
#plt.savefig('../figures/vc_amudiff_ext_15.png', dpi=500)
plt.show()


