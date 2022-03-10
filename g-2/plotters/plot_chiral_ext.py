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
from gvar import log
import sys

vc_amus_diff = vc_amu["ms"][2.5][1] 
vc_up_phys_diff = vc_phys15['up-charge'] - vc_phys15['up-nocharge']
vc_down_phys_diff = vc_phys15['down-charge'] - vc_phys15['down-nocharge']
ms_tuning = gv.load('../pseudoscalar/mistune_ps_fine.p')

#xdata = np.array( [3, 5, 7] )
xdata = np.array( [0.0362, 0.0364, 0.0366, 0.0368] )
#ydata = {"up": vc_amu_up, "down":vc_amu_down}
ydata = {   "Q=0":[y**2 for y in ms_tuning['E:n']], 
        "Q=0.101":[y**2 for y in ms_tuning['E:d']], 
        "Q=0.202":[y**2 for y in ms_tuning['E:u']]  }
print(ydata)
def fcn(x,p):
    out = {}
    for pt in x:
        out['Q=0']       = p['Q=0'][0] * (1 + p['Q=0'][1]*pt + p['Q=0'][2]*pt**2 )  
        out['Q=0.101']   = p['Q=0.101'][0] * (1 + p['Q=0.101'][1]*pt )  
        out['Q=0.202']   = p['Q=0.202'][0] * (1 + p['Q=0.202'][1]*pt )  
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
prior['Q=0'] = [gv.gvar('0(1)'),gv.gvar('0(1)'),gv.gvar(0,1)]#,gv.gvar(0,1)]
prior['Q=0.101'] = [gv.gvar('0(1)'),gv.gvar('0(1)')] #,gv.gvar(0,1),gv.gvar(0,1)]
prior['Q=0.202'] = [gv.gvar('0(1)'),gv.gvar('0(1)')] #,gv.gvar(0,1),gv.gvar(0,1)]
fit = lsq.nonlinear_fit(data=(xdata, ydata), prior=prior, fcn=fcn)
print("Fit results")
print(fit)

#ext_up = fcn(1, fit.p)['up']
#ext_down = fcn(1, fit.p)['down']
#ext_down_test = fcn(0.999294, fit.p)['down']
#gv.dump(extrapolated, 'amuu.p')
#print(ext_down)
#print(ext_down_test)

sys.exit(0)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

#fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
fig, ax1 = plt.subplots()

fitrange = np.arange(0.5,7.1,0.01)
fitline = [fcn(x,fit.p) for x in fitrange]

## plot squared PS masses


#plt.errorbar([x+0.00001 for x in xdata], [y.mean for y in Mps_Q202], xerr=0, yerr=[y.sdev for y in Mps_Q202], fmt='h', color='red', label='$Q=2e/3$', marker='x')
#plt.errorbar([x-0.00001 for x in xdata], [y.mean for y in Mps_Q101], xerr=0, yerr=[y.sdev for y in Mps_Q101], fmt='h', color='green', label='$Q=e/3$', marker='+')
#plt.errorbar(xdata, [y.mean for y in Mps_Q0], xerr=0, yerr=[y.sdev for y in Mps_Q0], fmt='h', color='black', label='$Q=0$')

#ax1.errorbar([x for x in xdata] ,[y.mean for y in ydata["up"]],xerr=0,yerr=[y.sdev for y in ydata["up"]] ,fmt='h',color='r',lw=1, label="$Q=2e/3$")
ax1.errorbar([x for x in xdata] ,[y.mean for y in ydata["down"]],xerr=0,yerr=[y.sdev for y in ydata["down"]] ,fmt='h',color='blue',lw=1, label="$Q=e/3$")
ax1.errorbar(27 , vc_amus_diff.mean ,xerr=0,yerr=vc_amus_diff.sdev ,fmt='h',color='blue',lw=1, label="$m_s$")

# extrapolated 
#ax1.errorbar(2/3+0.01, ext_up.mean,xerr=0,yerr=ext_up.sdev,fmt='h',mfc='none', label='$m_u$ linear extrapolation', color='red', lw=1, marker='*')
ax1.errorbar(4/3+0.01, ext_down.mean,xerr=0,yerr=ext_down.sdev,fmt='h',mfc='none', label='$m_d$ linear extrapolation', color='blue', lw=1, marker='*')

#ax1.plot([p for p in fitrange],[p['up'].mean for p in fitline],'--',color='r')
#ax1.fill_between([p for p in fitrange],[p['up'].mean-p['up'].sdev for p in fitline],[p['up'].mean+p['up'].sdev for p in fitline],alpha=0.5,lw=0,color='r')
#ax1.plot([p for p in fitrange],[p['down'].mean for p in fitline],'--',color='blue')
ax1.fill_between([p for p in fitrange],[p['down'].mean-p['down'].sdev for p in fitline],[p['down'].mean+p['down'].sdev for p in fitline],alpha=0.5,lw=0,color='blue')

## phys values
#ax1.errorbar(2/3-0.02, vc_up_phys_diff.mean,xerr=0,yerr=vc_up_phys_diff.sdev,fmt='h',label='$m_u, m_d$ sim $t^* = 1.5$fm', color='black', lw=1, marker='s', markersize=4)
#ax1.errorbar(4/3-0.02, vc_down_phys_diff.mean,xerr=0,yerr=vc_down_phys_diff.sdev,fmt='h',color='black', lw=1, marker='s', markersize=4)

ax1.legend(frameon=False, loc='lower right')
#ax1.legend(frameon=False, loc='upper right')
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=2/3, color='blue', linestyle='--', lw=0.5)
#ax1.axvline(x=4/3, color='blue', linestyle='--', lw=0.5)

ax1.set_xticks([2/3,4/3,3,5,7])
ax1.set_xticklabels(['$m_u$','$m_d$','3','5','7'])
#ax1.set_ylabel(r'$\delta a^{\tiny\mbox{HVP}}_{\mu}$', rotation=0, fontsize=15, labelpad=10)
ax1.set_xlabel(r'$\frac{m_q}{m_l}$', fontsize=25, labelpad=0)
#ax2.set_xlabel(r'$am_q$', fontsize=25, labelpad=20)
#plt.ylabel(r'$\frac{a_{\mu}^{\mathrm{qcd+qed}}}{a_{\mu}^{\mathrm{qcd}}}$', rotation=0, labelpad=25, fontsize=20)
#plt.ylabel(r'$a^2M^2_{\eta}$', rotation=0, labelpad=25, fontsize=20)
ax1.set_title('0.15fm ensemble - 1844 configs')
#fig.text(0.00, 0.55, r'$\delta a^{\small\mbox{HVP}}_{\mu}$', va='center', fontsize=24)
#ax1.set_ylim(bottom=-3e-10, top=0.8e-10)

plt.tight_layout()
#plt.savefig('../figures/vc_amudiff_ext_15.png', dpi=500)
plt.show()

plt.close()

