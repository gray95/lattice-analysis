import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt

mphys    = 0.002426				# isospin averaged
seven_ml = 0.01698/mphys 
five_ml =  0.01213/mphys
three_ml = 0.007278/mphys

#  	t*=10
amu_phys = gv.gvar('5.01(36)')
amu_thr  = gv.gvar('4.96(12)')
amu_fve  = gv.gvar('4.662(80)')
amu_svn  = gv.gvar('4.372(61)')

rt_svn = gv.gvar('0.99955(36)')
rt_fve = gv.gvar('1.00002(66)')
rt_thr = gv.gvar('1.0023(16)')
rphys = gv.gvar('1.003(15)')

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('500(200)'),gv.gvar(0,50),gv.gvar(0,1),gv.gvar(0,1)]


## extrapolation
fit = lsq.nonlinear_fit(data=([mphys*three_ml/hbarc,mphys*five_ml/hbarc,mphys*seven_ml/hbarc],[amu_thr, amu_fve, amu_svn]),prior=prior,fcn=extrap)
print("Fit results")
print(fit)

print("physical value: ", extrap([mphys/hbarc], fit.p))

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')


fitrange = np.arange(0,(seven_ml*mphys)**2,0.0001)
fitline = extrap([p/hbarc for p in fitrange],fit.p)

plt.errorbar((mphys*seven_ml)**2,amu_svn.mean,xerr=0,yerr=amu_svn.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar((mphys*five_ml)**2,amu_fve.mean,xerr=0,yerr=amu_fve.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar((mphys*three_ml)**2,amu_thr.mean,xerr=0,yerr=amu_thr.sdev,fmt='h',mfc='none',color='r',lw=1)

# physical point
plt.errorbar(mphys**2,amu_phys.mean,xerr=0,yerr=amu_phys.sdev,fmt='h',mfc='none', label='physical point', color='b', lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower left')

plt.xlabel('$m_q^2$')
plt.ylabel(r'$a_{\mu}^{\mbox{QCD+QED}}$', rotation=0, labelpad=30)

#plt.ylim(bottom=.98)
#plt.xlim(left=0, right=10)

plt.tight_layout()
#plt.savefig('amu_ratio_chiExt.png', dpi=600)
plt.show()

plt.close()

