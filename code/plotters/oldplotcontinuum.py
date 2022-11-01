import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from results import onemp_m, pp_m, onemp_am, pp_am
from commonweal import ainv, hbarc
from parameters import save_res_name
import collections

onempm = onemp_m

a = [ 1/ainv[x] for x in onemp_m ]
aa = [ 1/ainv[x] for x in pp_m ]
onemp_m = [ onemp_m[y] for y in onemp_m ]
pp_m = [ pp_m[y] for y in pp_m ]

# correlation
#vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1]])

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 + p['a'][2]*(point)**4 ))#p['a'][2]*mistunings[i]))
    return res


prior = {}
prior['a'] = [gv.gvar('4(2)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]
mistunings = [0, 0, 3.4, 8.8, 0]

## a**2 extrapolation
fit_h = lsq.nonlinear_fit(data=([x.mean/hbarc for x in a], onemp_m),prior=prior,fcn=extrap)
fit_pp = lsq.nonlinear_fit(data=([x.mean/hbarc for x in aa], pp_m),prior=prior,fcn=extrap)
print("Fit results")
print(fit_h, fit_pp)

continuum_limit_h = extrap( [0], fit_h.p )
continuum_limit_pp = extrap( [0], fit_pp.p )
print(continuum_limit_h, continuum_limit_pp)

## save ext values
cont_masses = {'onemp':continuum_limit_h[0], 'pp':continuum_limit_pp[0]}

gv.dump(cont_masses, '../out/FINAL_RES', add_dependencies=True)


# error budget
#print("Correlators - %s\nFit model - %s"%(corrpath,save_fit))
#fit = pickle.load(bz2.BZ2File(save_fit, 'rb'))
#print("reduced chi-squared: ", fit.chi2/fit.dof)
#print("t0 is {}\nfit-range [{},{}]".format(t0, TFIT[0], TFIT[-1]))
load_res = gv.load(save_res_name)
print(load_res['fit.p'])
outputs = collections.OrderedDict()
outputs['cont onemp'] = continuum_limit_h[0]
outputs['cont pp'] = continuum_limit_pp[0]
inputs = collections.OrderedDict()
inputs['onemp'] = onemp_am
inputs['latt space'] = ainv
inputs['prior'] = prior
print(gv.fmt_values(outputs))
print(gv.fmt_errorbudget(outputs, inputs, colwidth=18, verify=True))

## Plot extrapolation 

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(-0.005,a[0].mean,0.0001)
fitline_h = extrap([p/hbarc for p in fitrange],fit_h.p)
fitline_pp = extrap([p/hbarc for p in fitrange],fit_pp.p)

plt.errorbar([i.mean**2 for i in a], [x.mean for x in onemp_m],xerr=0,yerr=[y.sdev for y in onemp_m],fmt='o',mfc='none',color='r',label=r'NO - $1^{-+}$ hybrid',lw=1)
plt.errorbar([i.mean**2 for i in aa], [x.mean for x in pp_m],xerr=0,yerr=[y.sdev for y in pp_m],fmt='x',mfc='none',color='g',label='O - $1^{++}$ parity partner',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline_h],'--',color='r')
plt.plot([p**2 for p in fitrange],[p.mean for p in fitline_pp],'--',color='g')

plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_h],[p.mean+p.sdev for p in fitline_h],alpha=0.3,lw=0,color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline_pp],[p.mean+p.sdev for p in fitline_pp],alpha=0.3,lw=0,color='g')

# physical point
plt.errorbar(0.0, continuum_limit_h[0].mean,xerr=0,yerr=continuum_limit_h[0].sdev,fmt='+',mfc='none', color='red', lw=1)
plt.errorbar(0.0, continuum_limit_pp[0].mean,xerr=0,yerr=continuum_limit_pp[0].sdev,fmt='+',mfc='none', color='green', lw=1)

#plt.text(0.003, 2.65, 'PRELIMINARY', fontsize=6, bbox={'facecolor': 'pink', 'alpha': 0.7, 'pad': 4})
plt.errorbar(0.0004, 3.51067,xerr=0.0,yerr=0.05/1000,fmt='*', label=r"Expt. $\chi_{c1}(1P)$", color='k', lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$', fontsize=15)
plt.ylabel('mass\n [GeV]', rotation=0, labelpad=30, fontsize=14)
#plt.title(r'Continuum extrapolation of $\bar{c}c$ mass')
plt.ylim(top=5, bottom=2.5)
plt.xlim(left=-0.005)


plt.tight_layout()
#plt.savefig('../../figures/TESThybrid_continuum_ext.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

