import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from amu_results import bmw_diff, giu_diff


a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
#a_fine = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))

a = [a_vcoarse, a_coarse]

# data hard wired
#amu       = []
amu_diff = [gv.gvar('-0.00060(86)'), gv.gvar('-0.0021(13)')]	# [vc,c]

# correlation
#vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1]])

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('0(1)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]

print(1/(a_vcoarse.mean/hbarc),1/(a_coarse.mean/hbarc))

## a**2 extrapolation
fit = lsq.nonlinear_fit(data=([a_vcoarse.mean/hbarc,a_coarse.mean/hbarc], amu_diff),prior=prior,fcn=extrap)
print("Fit results")
print(fit)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,a_vcoarse.mean,0.0001)
fitline = extrap([p/hbarc for p in fitrange],fit.p)

plt.errorbar([i.mean**2 for i in a], [x.mean for x in amu_diff],xerr=0,yerr=[y.sdev for y in amu_diff],fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.3,lw=0,color='r')


plt.errorbar(0,bmw_diff.mean,xerr=0,yerr=bmw_diff.sdev,fmt='s',color='b',label='BMW',lw=1)
plt.errorbar(0,giu_diff.mean,xerr=0,yerr=giu_diff.sdev,fmt='s',color='g',label='Giusti et al.',lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=10,loc='upper right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.ylabel(r'$\frac{\delta a_{\mu}}{a_{\mu}}$', rotation=0, labelpad=10, fontsize=20)

#plt.ylim(top=0.00)
plt.xlim(left=-0.0005)

plt.tight_layout()
plt.savefig('../figures/amu_qedcorrection_summary.png', dpi=500, bbox_inches="tight")
#plt.show()

plt.close()

