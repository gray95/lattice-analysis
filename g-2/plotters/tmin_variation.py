
import sys
sys.path.append('..')
import gvar as gv
import matplotlib.pyplot as plt

tmin = range(1,6)

c = {}
c[1] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin1.p'
c[2] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin2.p'
c[3] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin3.p'
c[4] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin4.p'
c[5] = '/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin5.p'

c_dd = gv.BufferDict()

for t in tmin:
   tmp = gv.load(c[t])
   c_dd[t] = [ (y2-y1) for (y2,y1) in zip(tmp['down-charge'],tmp['down-nocharge']) ]

print(c_dd)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':18})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)


plt.errorbar([t-0.2 for t in tmin], [c_dd[t][0].mean for t in tmin], yerr=[c_dd[t][0].sdev for t in tmin], color='black', marker='o', markersize=6, label=r'$\Delta Q=1/3, 3m_l$', fmt='h', alpha=0.5)
plt.errorbar([t-0.1 for t in tmin], [c_dd[t][1].mean for t in tmin], yerr=[c_dd[t][1].sdev for t in tmin], color='purple', marker='o', markersize=6, label=r'$\Delta Q=1/3, 5m_l$', fmt='h', alpha=0.5)
plt.errorbar([t for t in tmin], [c_dd[t][2].mean for t in tmin], yerr=[c_dd[t][2].sdev for t in tmin], color='blue', marker='o', markersize=6, label=r'$\Delta Q=1/3, 7m_l$', fmt='h', alpha=0.5)
plt.errorbar([t+0.1 for t in tmin], [c_dd[t][3].mean for t in tmin], yerr=[c_dd[t][3].sdev for t in tmin], color='red', marker='o', markersize=6, label=r'$\Delta Q=1/3, m_s$', fmt='h', alpha=0.5)


#plt.xticks(tmin, [r'$B{=}1$',r'$B{=}2$',r'$B{=}4$',r'$B{=8}$'])

#plt.xlabel('$m_q/m_l$', fontsize=15)
plt.ylabel(r'$\Delta a_{\mu}$', rotation=0, labelpad=25, fontsize=24)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('../figures/amudowndiff_c_tminvar.png', dpi=500, bbox_inches="tight")
plt.show()

