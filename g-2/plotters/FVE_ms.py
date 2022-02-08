
import matplotlib.pyplot as plt
from results import coarseFV_amus_diff, coarseFV_amus_rt, coarseFV_amus
from results import coarseFV2Q_amus_diff, coarseFV2Q_amus_rt
from results import coarseFV_phi_m_diff, coarseFV_phi_f_diff, coarseFV_phi_m_rt
from results import coarseFV_phi_f_rt, coarseFV2Q_phi_f_rt

import gvar as gv


fig, (ax1, ax2) = plt.subplots( 2, 1, sharex=True )
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':14})
plt.rc('axes', linewidth=0.5)
plt.rc('text', usetex=True)


L = [24, 40, 48]

y1data = coarseFV_phi_m_rt
y1mean = 1/3 * (y1data[0]+y1data[1]+y1data[2])

y2data = coarseFV_amus_rt
y2mean = 1/3 * (y2data[0]+y2data[1]+y2data[2])

y3data = coarseFV2Q_amus_rt
y3mean = 1/3 * (y3data[0]+y3data[1]+y3data[2])


ax1.errorbar( [1/x for x in L], [y.mean for y in y1data], xerr=0, yerr=[y.sdev for y in y1data], fmt='h', mfc='none', color='blue', label='this work' )
ax1.axhline(y=y1mean.mean, color='red', linestyle='--')
#ax1.set_ylabel(r'$\delta M_{\phi}$', rotation=0, labelpad=20, fontsize=20)
ax1.set_ylabel(r'$R^0_{\rm{QED}}[M_{\phi}]$', rotation=0, labelpad=40, fontsize=18)
#ax1.set_ylabel(r'$a^{(s)}_{\mu}$', rotation=0, labelpad=20, fontsize=18)


ax2.errorbar( [1/x for x in L], [y.mean for y in y2data], xerr=0, yerr=[y.sdev for y in y2data], fmt='h', mfc='none', color='blue', label=r'$Q=1/3$' )
ax2.errorbar( [1/x for x in L], [y.mean for y in y3data], xerr=0, yerr=[y.sdev for y in y3data], fmt='h', mfc='none', color='green', label=r'$Q=2/3$' )
ax2.axhline(y=y2mean.mean, color='red', linestyle='--')
ax2.axhline(y=y3mean.mean, color='orange', linestyle='--')
#ax2.set_ylabel(r'$\delta f_{\phi}$', rotation=0, labelpad=20, fontsize=20)
ax2.set_ylabel(r'$R^0_{\rm{QED}}[a^{(s)}_{\mu}]$', rotation=0, labelpad=40, fontsize=18)
#ax2.set_ylabel(r'$R^0_{\rm{QED}}[f_{\phi}]$', rotation=0, labelpad=40, fontsize=18)
#ax2.set_ylabel(r'$\delta a^{(s)}_{\mu}$', rotation=0, labelpad=20, fontsize=18)
ax2.set_xlabel(r'$1/L_s$', fontsize=20)

#ax1.legend(frameon=False, loc='upper left')
ax2.legend(frameon=False, loc='center')

plt.tight_layout()

plt.savefig('../figures/update_FVE_ms_phiM.png', dpi=500)

plt.show()



