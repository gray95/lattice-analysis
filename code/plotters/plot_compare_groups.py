# script to plot various quantities 

import matplotlib.pyplot as plt
import gvar as gv
from commonweal import a, onemp_a, onemp_mass, plym_onemp_gev

fig, ax = plt.subplots(figsize=(13,7))
ax.set_yscale('linear')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

ax.errorbar([x.mean**2 for x in a], [y.mean for y in plym_onemp_gev], yerr=[y.sdev for y in plym_onemp_gev], xerr=0, fmt='h', mfc='none', color='red', label='this work')

ax.errorbar(onemp_a['hadspec']**2, onemp_mass['hadspec'].mean, xerr=0, yerr=onemp_mass['hadspec'].sdev, fmt='h',mfc='none',color='blue',lw=3, label='HadSpec [1610.01073]', alpha=0.8)

ax.errorbar(onemp_a['ma']**2, onemp_mass['ma'].mean, xerr=0, yerr=onemp_mass['ma'].sdev, fmt='h',mfc='none',color='black',lw=3, label='Ma et al. [1910.09819]', alpha=0.8)

ax.errorbar(onemp_a['regensburg']**2, onemp_mass['regensburg'].mean, xerr=0, yerr=onemp_mass['regensburg'].sdev, fmt='h',mfc='none',color='green',lw=3, label='Bali et al. [1110.2381]', alpha=0.8)

# DECAY THRESHOLDS
plt.plot([0.0 , 0.025 ] , [ 2*1.865 , 2*1.865  ] , "--", color = "black" )
plt.text(0.001 , 2*1.865 +0.05  , r"D$\overline{D}$ expt. threshold")

# two P-wave 
plt.text(0.001 , 2*2.343 +0.05  , r"$D_0^{\star}(2300)\overline{D_0}^{\star}(2300)$ expt. threshold")
plt.plot([0.0 , 0.025 ] , [ 2*2.343 , 2*2.343  ] , "--", color = "black" )
##------------------------------------------------------------------------------

ax.set_xlabel('$a^2$ [fm]', fontsize=14, labelpad=0)
ax.set_ylabel('$1^{-+}$ mass\n [GeV]', fontsize=14, rotation=0, labelpad=50)
ax.legend(frameon=False,fontsize=15,loc='lower right')
plt.ylim(top=5, bottom=4.0)
plt.xlim(left=0.0)
#plt.savefig('../figures/onemp_summary_plot.png', dpi=500, bbox_inches="tight")
#plt.tight_layout()
plt.show()

plt.close()
