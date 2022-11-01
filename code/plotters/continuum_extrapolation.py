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

a = [1/ainv[x] for x in onemp_m ]
onemp_m = [onemp_m[y] for y in onemp_m ]
pp_m = [pp_m[y] for y in pp_m ]

xdata = {'a':a}
ydata = {'onemp':onemp_m, 'pp':pp_m}

print(xdata)
print(ydata)

mistunings = {'vc':0, 'c':0, 'f':0, 'c-10':3.4, 'f-5':8.8} 
m = [mistunings[i] for i in onemp_am]

def extrap(p):
    res = {'onemp':[], 'pp':[]}
    n0, n1, n2, n3 = p['NO']
    o0, o1, o2, o3 = p['O']
    aa = p['a']
    Y = 1     # lambda QCD
#    o2 = 0.0
    for (a,dm) in zip(aa,m):
        res['onemp'].append(n0*(1 + n1*(Y*a)**2 + n2*(Y*a)**4 + n3*dm ))
        res['pp'].append(   o0*(1 + o1*(Y*a)**2 + o2*(Y*a)**4 + o3*dm ))
    return res

def make_prior(x):
    prior = gv.BufferDict()
    prior['NO'] =  [gv.gvar('4(2)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]
    prior['O'] =   [gv.gvar('4(2)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]
    prior['a'] = x['a']
    return prior


prior = make_prior(xdata)
fit = lsq.nonlinear_fit(data=ydata, prior=prior, fcn=extrap)
print("fit results")
print(fit)


continuum = [gv.gvar('0(0)')]
X = make_prior({'a':continuum})
X['NO'] = fit.p['NO']
X['O']  = fit.p['O']

ext = extrap(X)
print(ext)
## error analysis
#print("Correlators - %s\nFit model - %s"%(corrpath,save_fit))
#fit = pickle.load(bz2.BZ2File(save_fit, 'rb'))
#print("reduced chi-squared: ", fit.chi2/fit.dof)
#print("t0 is {}\nfit-range [{},{}]".format(t0, TFIT[0], TFIT[-1]))
load_res = gv.load(save_res_name)
print(load_res['fit.p'])
outputs = collections.OrderedDict()
outputs['cont onemp'] = ext['onemp'][0]
outputs['cont pp'] = ext['pp'][0]
inputs = collections.OrderedDict()
inputs['onemp'] = onemp_am
inputs['latt space'] = a
inputs['a-->0'] = {'NO':prior['NO'], 'O':prior['O']}
print(gv.fmt_values(outputs))
print(gv.fmt_errorbudget(outputs, inputs, colwidth=18, verify=True))

sys.exit(0)

## Plot extrapolation 

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(-0.005,a[0].mean,0.001)
m = [0 for i in np.arange(-0.005,a[0].mean,0.001)]
X = make_prior({'a':fitrange})
X['NO'] = fit.p['NO']
X['O']  = fit.p['O']
fitline = extrap(X)

plt.errorbar([i.mean**2 for i in a], [x.mean for x in onemp_m],xerr=0,yerr=[y.sdev for y in onemp_m],fmt='o',mfc='none',color='r',label=r'NO - $1^{-+}$ hybrid',lw=1)
plt.errorbar([i.mean**2 for i in a], [x.mean for x in pp_m],xerr=0,yerr=[y.sdev for y in pp_m],fmt='x',mfc='none',color='g',label='O - $1^{++}$ parity partner',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline['onemp']],'--',color='r')
plt.plot([p**2 for p in fitrange],[p.mean for p in fitline['pp']],'--',color='g')

plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline['onemp']],[p.mean+p.sdev for p in fitline['onemp']],alpha=0.3,lw=0,color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline['pp']],[p.mean+p.sdev for p in fitline['pp']],alpha=0.3,lw=0,color='g')

# physical point
plt.errorbar(0.0, ext['onemp'][0].mean,xerr=0,yerr=ext['onemp'][0].sdev,fmt='+',mfc='none', color='red', lw=1)
plt.errorbar(0.0, ext['pp'][0].mean,xerr=0,yerr=ext['pp'][0].sdev,fmt='+',mfc='none', color='green', lw=1)

plt.errorbar(0.0004, 3.51067,xerr=0.0,yerr=0.05/1000,fmt='*', label=r"Expt. $\chi_{c1}(1P)$", color='k', lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower right')

plt.xlabel('$a^2\ [\mbox{GeV}^{-2}]$', fontsize=15)
plt.ylabel('mass\n[GeV]', rotation=0, labelpad=15, fontsize=14)
#plt.title(r'Continuum extrapolation of $\bar{c}c$ mass')
plt.ylim(top=5, bottom=2.5)
plt.xlim(left=-0.005)


plt.tight_layout()
plt.savefig('../../figures/onemp_continuum_ext.png', dpi=500, bbox_inches="tight")
plt.show()

plt.close()

