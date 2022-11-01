
## visualise all the QED corrections in a compact way

import gvar as gv
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from results import vc_down_diff, c_down_diff, f_down_diff
from results import vc_up_diff, c_up_diff, f_up_diff

print(vc_up_diff)
print(f_up_diff)
sys.exit(0)
masses = [3, 5, 7, 27.8]

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':18})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

plt.errorbar(range(4), [y.mean for y in vc_down_diff], yerr=[y.sdev for y in vc_down_diff], fmt='h', color='blue', marker='o', markersize=6, label=r'$\Delta Q=1/3,0.15$fm', mfc='none')
plt.errorbar([i-0.1 for i in range(4)], [y.mean for y in c_down_diff], yerr=[y.sdev for y in c_down_diff], fmt='h', color='blue', marker='s', markersize=6, label=r'$\Delta Q=1/3,0.12$fm', mfc='none')
plt.errorbar([i-0.2 for i in range(4)], [y.mean for y in f_down_diff], yerr=[y.sdev for y in f_down_diff], fmt='h', color='blue', marker='D', markersize=6, label=r'$\Delta Q=1/3,0.09$fm', mfc='none')

plt.errorbar(range(4), [y.mean for y in vc_up_diff], yerr=[y.sdev for y in vc_up_diff], fmt='h', color='red', marker='o', markersize=6, label=r'$\Delta Q=2/3,0.15$fm', alpha=0.5)
plt.errorbar([i-0.1 for i in range(4)], [y.mean for y in c_up_diff], yerr=[y.sdev for y in c_up_diff], fmt='h', color='red', marker='s', markersize=6, label=r'$\Delta Q=2/3,0.12$fm', alpha=0.5)
plt.errorbar([i-0.2 for i in range(4)], [y.mean for y in f_up_diff], yerr=[y.sdev for y in f_up_diff], fmt='h', color='red', marker='D', markersize=6, label=r'$\Delta Q=2/3,0.09$fm', alpha=0.5)

plt.axhline(y=0, linestyle='--', linewidth=0.5)


plt.xticks(range(4), [r'$3m_l$',r'$5m_l$',r'$7m_l$',r'$m_s$'])

#plt.xlabel('$m_q/m_l$', fontsize=15)
plt.ylabel(r'$\Delta a_{\mu}$', rotation=0, labelpad=15, fontsize=24)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('../figures/amudiffdata_summary.png', dpi=500, bbox_inches="tight")
#plt.show()

