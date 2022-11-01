# plot the conventional meson properties - mass and decay constant
# qed corrections
# eta mass and decay constant plus corrections
# add in Bijnens number 

import sys
sys.path.append('..')
import matplotlib.pyplot as plt
import gvar as gv
from commonweal import hbarc, w0, w0overa

a = [ w0/(hbarc*w0overa[s]) for s in ['vc','c','f'] ]

ams = [ 0.0677, 0.0527, 0.0364 ]    # strange quark in lattice units

vc = gv.load('../pseudoscalar/fits/ps_vcoarse.p')
c = gv.load('../pseudoscalar/fits/ps_coarse.p')
f = gv.load('../pseudoscalar/fits/ps_fine.p')

m_etas = [ 1000*X['E:n'][-1]/a for (X,a) in zip([vc, c, f],a) ]
f_etas = [ 1000*(8**0.5)*X['a:n'][-1]*(1/a)*ms*(X['E:n'][-1]**(-3/2)) for (X,ms,a) in zip([vc,c,f],ams,a) ] 

m_ratio = [ X['E:d'][-1]/X['E:n'][-1] for X in [vc,c,f] ]
f_ratio = [ (X['a:d'][-1]/X['a:n'][-1])*(X['E:d'][-1]/X['E:n'][-1])**(-1.5) for X in [vc,c,f] ]

m_diff = [ 1000*(X['E:d'][-1]-X['E:n'][-1])/a for (X,a) in zip([vc, c, f],a) ]
f_diff = [ 1000*(8**0.5)*1/a*ms*(X['a:d'][-1]*X['E:d'][-1]**(-3/2)-X['a:n'][-1]*X['E:n'][-1]**(-3/2)) for (X,ms,a) in zip([vc,c,f],ams,a) ] 

print(m_diff)
print(f_diff)
print(m_etas)
print(f_etas)
xdata = [(hbarc*x)**2 for x in a]  # in fm^2

#### PLOTTING BEGINS ####

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.errorbar( [x.mean for x in xdata], [y.mean for y in f_diff], xerr=[x.sdev for x in xdata], yerr=[y.sdev for y in f_diff], fmt='h', color='red',lw=1, label='this work', alpha=0.8)
#ax1.errorbar( 0, 689.89, yerr=0.49, color='blue', fmt='h', label='BMW (2020)', alpha=0.7)
ax2.errorbar( [x.mean for x in xdata], [y.mean for y in m_diff], xerr=[x.sdev for x in xdata], yerr=[y.sdev for y in m_diff], fmt='h', color='red',lw=1, label='this work', alpha=0.8)


ax1.set_xlim(left=-0.001)
ax2.set_xlabel('$a^2$[fm]', fontsize=25, labelpad=20)
ax1.set_ylabel('$\Delta f_{\eta_s}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax1.yaxis.set_label_coords(-0.2,0.25)
ax2.yaxis.set_label_coords(-0.2,0.25)
ax2.set_ylabel('$\Delta M_{\eta_s}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax1.legend(frameon=False,fontsize=12,loc='lower left')
#ax1.legend(frameon=False,fontsize=12,loc='upper left')

plt.tight_layout()
plt.savefig('../figures/etas_qedproperties.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()
