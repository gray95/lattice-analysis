# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv

## COARSE
# qcd+qed values
c_amu = [gv.gvar('5.31(13)'), gv.gvar('4.872(77)'), gv.gvar('4.549(58)')]
# rt of qcdqed/qcd
c_rt = [gv.gvar('0.9980(12)'), gv.gvar('0.99845(54)'), gv.gvar('0.99862(31)')]
# diff - qcdqed-qcd/qcd
c_diff = [gv.gvar('-0.0020(12)'), gv.gvar('-0.00155(54)'), gv.gvar('-0.00138(31)')]

## VERY COARSE
vc_amu = [gv.gvar('5.09(10)'), gv.gvar('4.727(68)'), gv.gvar('4.429(55)')]
vc_rt  = [gv.gvar('0.9997(12)'), gv.gvar('0.99895(51)'), gv.gvar('1.0009(79)')]
vc_diff= [gv.gvar('-0.0003(12)'), gv.gvar('-0.00105(51)'), gv.gvar('0.0009(79)')]

vc_amu_phys = gv.gvar('5.49(18)')		# combined a and b stream
vc_rt_phys  = gv.gvar('1.0057(65)')
vc_diff_phys= gv.gvar('0.0057(65)')

## BMW numbers
bmw_amuqcd = gv.gvar('633.7(2.1)')
bmw_amuqcdqed = bmw_amuqcd + gv.gvar('-1.23(40)') + gv.gvar('6.60(63)')
bmw_rt = bmw_amuqcdqed/bmw_amuqcd
bmw_diff = (bmw_amuqcdqed-bmw_amuqcd)/bmw_amuqcd

###############################################################################
## Plotting begins here##

fig, ax = plt.subplots(figsize=(10,7))
ax.set_yscale('linear')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

#ax.errorbar(1.05, (bmw_amuqcdqed.mean/100), xerr=0, yerr=(bmw_amuqcdqed.sdev)/100, fmt='h',mfc='none',color='black',lw=2, label='BMW - strong + qed isospin breaking included', alpha=0.8)
ax.errorbar(1, vc_amu_phys.mean, xerr=0, yerr=vc_amu_phys.sdev, fmt='h',mfc='none',color='green',lw=1, label='very coarse physical point', alpha=0.8)

ax.errorbar([3, 5, 7], [x.mean for x in vc_amu], xerr=0, yerr=[x.sdev for x in vc_amu], fmt='h',mfc='none',color='red',lw=1, label='very coarse-185 configs', alpha=0.8)

ax.errorbar([3.05, 5.05, 7.05], [x.mean for x in c_amu], xerr=0, yerr=[x.sdev for x in c_amu], fmt='h',mfc='none',color='blue',lw=1.5, label='coarse-56 configs', alpha=0.8)

ax.set_xlabel(r'mass / ml', fontsize=18)
#ax.set_ylabel(r'$a_{\mu}$x$10^{-8}$', fontsize=14)
ax.set_ylabel('$a_{\mu}^{[qcd+qed]}$\n\nx$10^{-8}$', fontsize=24, rotation=0, labelpad=50)
ax.legend(frameon=False,fontsize=15,loc='upper right')
#plt.ylim(top=6)
#plt.xlim(left=1)
plt.savefig('../figures/g-2_amu.png', dpi=500, bbox_inches="tight")
plt.tight_layout()
plt.show()

plt.close()
