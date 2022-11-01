
import sys
sys.path.append('..')
import gvar as gv
import matplotlib.pyplot as plt

bins = [1, 2, 4, 8]

vc = {}
vc[1] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_vc_bin1_tmin2.p'
vc[2] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_vc_bin2_tmin2.p'
vc[4] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_vc_bin4_tmin2.p'
vc[8] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_vc_bin8_tmin2.p'
vc_dd = gv.BufferDict()

for B in bins:
   tmp = gv.load(vc[B])
   vc_dd[B] = [ (y2-y1) for (y2,y1) in zip(tmp['down-charge'],tmp['down-nocharge']) ]

print(vc_dd)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':18})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)


plt.errorbar([b-0.05 for b in bins], [vc_dd[B][0].mean for B in bins], yerr=[vc_dd[B][0].sdev for B in bins], color='black', marker='o', markersize=6, label=r'$\Delta Q=1/3, 3m_l$', fmt='h', alpha=0.5)
plt.errorbar([b for b in bins], [vc_dd[B][1].mean for B in bins], yerr=[vc_dd[B][1].sdev for B in bins], color='purple', marker='o', markersize=6, label=r'$\Delta Q=1/3, 5m_l$', fmt='h', alpha=0.5)
plt.errorbar([b+0.05 for b in bins], [vc_dd[B][2].mean for B in bins], yerr=[vc_dd[B][2].sdev for B in bins], color='blue', marker='o', markersize=6, label=r'$\Delta Q=1/3, 7m_l$', fmt='h', alpha=0.5)
plt.errorbar([b+0.1 for b in bins], [vc_dd[B][3].mean for B in bins], yerr=[vc_dd[B][3].sdev for B in bins], color='red', marker='o', markersize=6, label=r'$\Delta Q=1/3, m_s$', fmt='h', alpha=0.5)
plt.xticks(bins, [r'$B{=}1$',r'$B{=}2$',r'$B{=}4$',r'$B{=8}$'])

#plt.xlabel('$m_q/m_l$', fontsize=15)
plt.ylabel(r'$\Delta a_{\mu}$', rotation=0, labelpad=15, fontsize=24)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('../figures/amudowndiff_vc_binsize.png', dpi=500, bbox_inches="tight")
plt.show()

