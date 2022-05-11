# calculate and plot the conventional meson properties
# phi mass and decay constant
# qed corrections
# eta mass and decay constant plus corrections
## add in Bijnens number 

import matplotlib.pyplot as plt
import gvar as gv
from commonweal import hbarc, w0, w0overa, ZV, ZVqedd
import sys
import dill as pickle
import bz2
import lsqfit

a = [ w0/(hbarc*w0overa[s]) for s in ['vc','c','f'] ]
Zv = [ ZV[L] for L in ['vc','c','f'] ]
Zvd = [ ZV[L]*ZVqedd[L] for L in ['vc','c','f'] ]

vc = '../fits/fitobj/vc_FIT_bin2.pbz2'
c  = '../fits/fitobj/c_FIT_bin2.pbz2'
f  = '../fits/fitobj/f_FIT_bin1.pbz2'

dE = {'n':[], 'd':[]}
amp = {'n':[], 'd':[]}

m = 'ms'
#for L in [vc, c, f]:
#    fit = pickle.load(bz2.BZ2File(L, 'rb'))
#    dE['n'].append(fit.p[m+':dE:n'][0])
#    dE['d'].append(fit.p[m+':dE:d'][0])
#    amp['n'].append(fit.p[m+':a:n'][0])
#    amp['d'].append(fit.p[m+':a:d'][0])

#phi_mass = [ 1000*m/a for (m,a) in zip(dE['n'],a) ]
#phi_f    = [ 1000*Zv*ampl*((2/m)**0.5)*(1/a) for (Zv,ampl,m,a) in zip(Zv,amp['n'],dE['n'],a) ]

#phi_mdiff = [ 1000*(m1-m2)/a for (m1,m2,a) in zip(dE['d'],dE['n'],a) ]
#phi_fdiff = [ 1000*2**0.5*(Zvd*a1*(m1**-0.5) - Zv*a2*(m2**-0.5)) for (Zvd,Zv,a1,a2,m1,m2) in zip(Zvd,Zv,amp['d'],amp['n'],dE['d'],dE['n']) ]

#phi_mratio = [ m1/m2 for (m1,m2) in zip(dE['d'],dE['n']) ]
#phi_fratio = [ (Zvd/Zv)*(a1/a2)*((m1/m2)**(-0.5)) for (Zvd,Zv,a1,a2,m1,m2) in zip(Zvd,Zv,amp['d'],amp['n'],dE['d'],dE['n']) ]

phi_mass = gv.gvar(['1037.2(5.5)', '1031.5(5.6)', '1025.3(6.2)'])
phi_mdiff = gv.gvar(['0.246(89)', '0.3(1.1)', '0.3(3.5)'])
phi_mratio = gv.gvar(['1.000237(86)', '1.0003(10)', '1.0003(34)'])

phi_f    = gv.gvar(['241.9(1.5)', '237.8(2.1)', '236.7(3.9)'])
phi_fdiff = gv.gvar(['0.035(87)', '0.03(88)', '0.02(2.00)'])
phi_fratio = gv.gvar(['1.00019(47)', '1.0002(60)', '1.000(19)'])

xdata = [(hbarc*x)**2 for x in a]  # in fm^2

#### PLOTTING BEGINS ####

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)


ax1.errorbar( [x.mean for x in xdata], [y.mean for y in phi_mass], xerr=[x.sdev for x in xdata], yerr=[y.sdev for y in phi_mass], fmt='h' ,color='red',lw=1, label='this work', alpha=0.7)
ax1.errorbar( 0, 1019.461, yerr=0.016, fmt='h', color='black', label='exp', alpha=0.9 )
ax1.errorbar( 0, 1023.2, yerr=4.2, fmt='h', color='blue', label='HPQCD (2017)', alpha=0.9 )



ax2.errorbar( [x.mean for x in xdata], [y.mean for y in phi_f], xerr=[x.sdev for x in xdata], yerr=[y.sdev for y in phi_f], fmt='h' ,color='red',lw=1, label='this work', alpha=0.7)
ax2.errorbar( 0, 228.5, yerr=3.6, fmt='h', color='black', label='exp', alpha=0.9 )
ax2.errorbar( 0, 237.6, yerr=3.0, fmt='h', color='blue', label='HPQCD (2017)', alpha=0.9 )

ax1.set_xlim(left=-0.001)
ax2.set_xlabel('$a^2$[fm]', fontsize=25, labelpad=20)
ax1.set_ylabel('$M_{\phi}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax2.set_ylabel('$f_{\phi}$\n[MeV]', fontsize=20, rotation=0, labelpad=40)
ax1.yaxis.set_label_coords(-0.2,0.25)
ax2.yaxis.set_label_coords(-0.2,0.25)
#ax1.legend(frameon=False,fontsize=10,loc='upper left')
ax2.legend(frameon=False,fontsize=10,loc='lower right')

plt.tight_layout()
plt.savefig('../figures/phi_properties.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()
