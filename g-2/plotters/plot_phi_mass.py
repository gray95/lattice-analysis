# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from commonweal import vc_phi_mass, c_phi_mass, f_phi_mass
from commonweal import vc_fphi, c_fphi, f_fphi, hbarc


fig, ax = plt.subplots(figsize=(10,7))
ax.set_yscale('linear')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

phi_mass = [vc_phi_mass, c_phi_mass, f_phi_mass]
phi_decay_constant = [vc_fphi, c_fphi, f_fphi]

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))


a = [ a_vcoarse, a_coarse, a_fine]
print(a)

xdata = [x.mean**2 for x in a]

ydata = phi_decay_constant

ax.errorbar( xdata, [y.mean for y in ydata], xerr=0, yerr=[y.sdev for y in ydata], fmt='h',mfc='none',color='red',lw=1, label='this work', alpha=0.8)

ax.errorbar(0, 0.114, xerr=0, yerr=0, fmt='h',mfc='none',color='blue',lw=1.5, label='exp', alpha=0.8)

ax.set_xlabel(r'$a^2$', fontsize=25, labelpad=20)
#ax.set_ylabel(r'$\frac{\delta a_{\mu}^{(\mathrm{light})}}{a_{\mu}^{(\mathrm{light})}}$', fontsize=26, rotation=0, labelpad=40)
#plt.yticks(fontsize=12)#ax.set_yticklabels(fontsize=12)
ax.legend(frameon=False,fontsize=15,loc='upper left')
#plt.ylim(top=1.05, bottom = 1.01)
#plt.xlim(left=-0.05, right=1)
plt.savefig('../figures/phi_decayconstant.png', dpi=500, bbox_inches="tight")
plt.tight_layout()
plt.show()

plt.close()
