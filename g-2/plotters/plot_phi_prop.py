# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from results import vc_phi_m, c_phi_m, f_phi_m
from results import vc_phi_f, c_phi_f, f_phi_f
from commonweal import hbarc


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

phi_mass = [vc_phi_m, c_phi_m, f_phi_m]
phi_f = [vc_phi_f, c_phi_f, f_phi_f]

phi_mass_exp = gv.gvar('1.0195(5)')
phi_f_exp = gv.gvar('228.5(3.6)')

phi_mass_hpqcd = gv.gvar('1.0232(42)') 
phi_f_hpqcd = gv.gvar('237.6(2.9)')

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))


a = [ a_vcoarse, a_coarse, a_fine]

xdata = [x.mean**2 for x in a]
ydata = phi_mass
yydata = phi_f 

ax1.errorbar( xdata, [y.mean for y in ydata], xerr=0, yerr=[y.sdev for y in ydata], fmt='h',mfc='none',color='red',lw=1, label='this work', alpha=0.8)
ax1.errorbar(3e-4, phi_mass_hpqcd.mean, xerr=0, yerr=phi_mass_hpqcd.sdev, fmt='h',mfc='none',color='blue',lw=1.5, label='HPQCD [1703.05552]', alpha=0.8)
ax1.errorbar(-1e-4, phi_mass_exp.mean, xerr=0, yerr=phi_mass_exp.sdev, fmt='h',mfc='none',color='black',lw=1.5, label='exp', alpha=0.8)
#ax1.set_xticks([])

ax2.errorbar( xdata, [y.mean for y in yydata], xerr=0, yerr=[y.sdev for y in yydata], fmt='h',mfc='none',color='red',lw=1, label='this work', alpha=0.8)
ax2.errorbar(3e-4, phi_f_hpqcd.mean, xerr=0, yerr=phi_f_hpqcd.sdev, fmt='h',mfc='none',color='blue',lw=1.5, label='HPQCD [1703.05552]', alpha=0.8)
ax2.errorbar(-1e-4, phi_f_exp.mean, xerr=0, yerr=phi_f_exp.sdev, fmt='h',mfc='none',color='black',lw=1.5, label='exp', alpha=0.8)


ax2.set_xlabel('$a^2$[fm]', fontsize=25, labelpad=20)
ax1.set_ylabel('$M_{\phi}$\n[GeV]', fontsize=20, rotation=0, labelpad=40)
ax2.set_ylabel('$f_{\phi}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax1.legend(frameon=False,fontsize=12,loc='upper left')

plt.tight_layout()
#plt.savefig('../figures/phi_properties.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()
