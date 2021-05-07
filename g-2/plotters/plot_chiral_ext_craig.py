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
rt_phys = gv.gvar('1.003(15)')

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('1(1)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]


## extrapolation
##fit = lsq.nonlinear_fit(data=([mphys*three_ml/hbarc,mphys*five_ml/hbarc,mphys*seven_ml/hbarc],[amu_thr, amu_fve, amu_svn]),prior=prior,fcn=extrap)
fit = lsq.nonlinear_fit(data=([3,5,7],[rt_thr, rt_fve, rt_svn]),prior=prior,fcn=extrap)
print("Fit results")
print(fit)

print("physical value: ", extrap([1], fit.p))

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')


fitrange = np.arange(0,8,0.01)
fitline = extrap([p for p in fitrange],fit.p)

plt.errorbar(7,rt_svn.mean,xerr=0,yerr=rt_svn.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(5,rt_fve.mean,xerr=0,yerr=rt_fve.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(3,rt_thr.mean,xerr=0,yerr=rt_thr.sdev,fmt='h',mfc='none',color='r',lw=1)

# physical point
plt.errorbar(1,rt_phys.mean,xerr=0,yerr=rt_phys.sdev,fmt='h',mfc='none', label='physical point', color='b', lw=1)

plt.plot([p for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower left')

plt.xlabel('$m_q/m_l$')
plt.ylabel(r'$a_{\mu}^{\mbox{\tiny QCD+QED}}$', rotation=0, labelpad=35)

#plt.ylim(bottom=.98)
#plt.xlim(left=0, right=10)

plt.tight_layout()
plt.savefig('amu_ratio_chiExt.png', dpi=600)
plt.show()

plt.close()

