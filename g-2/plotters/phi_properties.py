# calculate and plot the conventional meson properties
# phi mass and decay constant
# qed corrections
# eta mass and decay constant plus corrections
## add in Bijnens number 

import matplotlib.pyplot as plt
import gvar as gv
from commonweal import hbarc
import sys


pickle23 = {"encoding":"latin1"}    # for compatibility btwn python2/3

vc = gv.load('../fits/fine/ms.p', **pickle23)

print(vc)

sys.exit(0)

p = fit.p
g = collections.OrderedDict()
for i in range(NMASSES):
    m = 'm'+str(i)
    for Q in ['n', 'u', 'd']:
        g[m+'E:'+Q] = np.cumsum(p[m+'dE:'+Q])[0]
        g[m+'a:'+Q] = p[m+'a:'+Q][0]
print(gv.tabulate(g, ncol=2))

phi_mass_diff = [vc_phi_m_diff, c_phi_m_diff]#, f_phi_m_diff]
phi_f_diff = [vc_phi_f_diff, c_phi_f_diff]#, f_phi_f_diff]

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse  = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
a_fine    = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))

a = [a_vcoarse, a_coarse, a_fine]

xdata = [x.mean**2 for x in a]
ydata = phi_mass_diff
yydata = phi_f_diff 


#### PLOTTING BEGINS ####

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)


ax1.errorbar( xdata, [y.mean for y in ydata], xerr=0, yerr=[y.sdev for y in ydata], fmt='h',mfc='none',color='red',lw=1, label='this work', alpha=0.8)
#ax1.set_xticks([])

ax2.errorbar( xdata, [y.mean for y in yydata], xerr=0, yerr=[y.sdev for y in yydata], fmt='h',mfc='none',color='red',lw=1, label='this work', alpha=0.8)
#ax2.errorbar(-1e-4, phi_f_exp.mean, xerr=0, yerr=phi_f_exp.sdev, fmt='h',mfc='none',color='black',lw=1.5, label='exp', alpha=0.8)

ax1.set_xlim(left=0)
ax2.set_xlabel('$a^2$[fm]', fontsize=25, labelpad=20)
#ax1.set_ylabel('$M_{\phi}$\n[GeV]', fontsize=20, rotation=0, labelpad=40)
#ax2.set_ylabel('$f_{\phi}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax1.legend(frameon=False,fontsize=12,loc='upper left')

plt.tight_layout()
#plt.savefig('../figures/phi_properties.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()
