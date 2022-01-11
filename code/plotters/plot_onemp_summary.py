import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from commonweal import plym_onemp_gev, pp_gev
from commonweal import a # list of lattice spacings coarsest --> finest

hbarc = 0.197326968


print(a)

# correlation
#vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1]])

def extrap(x,p):
    res = []
    lambda_gev = 0.5
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 / lambda_gev**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res


prior = {}
prior['a'] = [gv.gvar('4(1)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]


## a**2 extrapolation
fit_h = lsq.nonlinear_fit(data=([x.mean/hbarc for x in a], plym_onemp_gev),prior=prior,fcn=extrap)
fit_pp = lsq.nonlinear_fit(data=([x.mean/hbarc for x in a], pp_gev),prior=prior,fcn=extrap)
print("Fit results")
print(fit_h, fit_pp)

continuum_limit_h = extrap( [0], fit_h.p )
continuum_limit_pp = extrap( [0], fit_pp.p )
print(continuum_limit_h, continuum_limit_pp)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,a[0].mean,0.0001)
fitline_h = extrap([p/hbarc for p in fitrange],fit_h.p)
fitline_pp = extrap([p/hbarc for p in fitrange],fit_pp.p)

plt.errorbar([i.mean**2 for i in a], [x.mean for x in plym_onemp_gev],xerr=0,yerr=[y.sdev for y in plym_onemp_gev],fmt='h',mfc='none',color='r',label=r'This work',lw=1)
#plt.errorbar([i.mean**2 for i in a], [x.mean for x in pp_gev],xerr=0,yerr=[y.sdev for y in pp_gev],fmt='h',mfc='none',color='g',label='$1^{++}$ $\overline{c}c$ (parity partner)',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline_h],'--',color='r')
#plt.plot([p**2 for p in fitrange],[p.mean for p in fitline_pp],'--',color='g')

plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_h],[p.mean+p.sdev for p in fitline_h],alpha=0.3,lw=0,color='r')
#plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_pp],[p.mean+p.sdev for p in fitline_pp],alpha=0.3,lw=0,color='g')

# other groups
##plt.errorbar([0.12**2], [1.233+2.984],xerr=0,yerr=[16/1000],fmt='s',mfc='none',color='b',label=r'HadSpec 1204.5425, $m_\pi$=400 MeV',lw=1)
plt.errorbar([0.12**2], [1.233+2.984],xerr=0,yerr=[16/1000],fmt='o',mfc='none',color='b',lw=1)

plt.errorbar([0.12**2], [1.326+2.984],xerr=0,yerr=[23/1000],fmt='o',mfc='none',color='b',label=r'HadSpec 1610.01073, $m_\pi$=240,400 MeV',lw=1)

plt.errorbar([0.1145**2], [4.154],xerr=0,yerr=[54/1000],fmt='s',mfc='none',color='g',label=r'Bali et al. 1110.2381, $m_\pi$=1000 MeV, $n_f$=2 ',lw=1)

###  decay thresholds
plt.plot([0.0 , 0.025 ] , [ 2*1.865 , 2*1.865  ] , "--", color = "black" , lw=0.4)
plt.text(0.001 , 2*1.865 +0.05  , r"D$\overline{D}$", fontsize=6, alpha=0.5)

plt.plot([0.0 , 0.025 ] , [ (2.343+2.007), (2.343+2.007) ] , "--", color = "black", lw=0.4)
plt.text(0.001 , (2.343+2.007)+0.05  , r"$D_0^{\star}(2300)\overline{D^{\star}}$", fontsize=6, alpha=0.5)

plt.plot([0.0 , 0.025 ] , [ (2.422+1.865), (2.422+1.865) ] , "--", color = "black",lw=0.4)
plt.text(0.001 , (2.422+1.865)-0.15  , r"$D_1\overline{D}$", fontsize=6, alpha=0.5)

# two P-wave 
plt.text(0.001 , 2*2.343 +0.05  , r"$D_0^{\star}(2300)\overline{D_0}^{\star}(2300)$", fontsize=6, alpha=0.5)
plt.plot([0.0 , 0.025 ] , [ 2*2.343 , 2*2.343  ] , "--", color = "black", lw=0.4 )
#---------------------------------------------------------------------------

plt.text(0.016, 3.9, 'PRELIMINARY', fontsize=6, bbox={'facecolor': 'pink', 'alpha': 0.7, 'pad': 4})

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=10)
plt.ylabel('mass\n [GeV]', rotation=0, labelpad=10, fontsize=10)
plt.title(r'Summary of $\overline{c}c$ $1^{-+}$ masses')
plt.ylim(top=5, bottom=3.3)
plt.xlim(left=-0.0002)


plt.tight_layout()
plt.savefig('../../figures/onemp_summary.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

