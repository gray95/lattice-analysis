import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt

# vc t*=10
vc_mphys    = 0.002426				# isospin averaged
vc_seven_ml = 0.01698 
vc_five_ml =  0.01213
vc_three_ml = 0.007278
vc_rt_7ml = gv.gvar('0.99955(36)') #	old single pol data
vc_rt_5ml = gv.gvar('0.99895(51)')
vc_rt_3ml = gv.gvar('0.9997(12)')
vc_rphys = gv.gvar('1.0057(65)')	# combined a and b stream
vc_diff_3ml = gv.gvar('-0.0003(12)') 
vc_diff_5ml = gv.gvar('-0.00105(51)')
vc_diff_7ml = gv.gvar('-0.00045(36)')
vc_diff_phys= gv.gvar('0.0057(65)')

# coarse t*=13
c_ml      = 0.00184
c_3ml     = 0.00552
c_5ml     = 0.0092
c_7ml     = 0.01288
c_rt_3ml = gv.gvar('0.9980(12)')
c_rt_5ml = gv.gvar('0.99845(54)')
c_rt_7ml = gv.gvar('0.99862(31)')
c_diff_3ml = gv.gvar('-0.0020(12)')
c_diff_5ml = gv.gvar('-0.00155(54)')
c_diff_7ml = gv.gvar('-0.00138(31)')

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point) )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('0(1)'),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]


## extrapolation
fit = lsq.nonlinear_fit(data=([3,5,7],[c_diff_3ml, c_diff_5ml, c_diff_7ml]),prior=prior,fcn=extrap)
print("Fit results")
print(fit)
extrapolated = extrap([1], fit.p)
print("physical value: ", extrapolated)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')


fitrange = np.arange(0,8,0.01)
fitline = extrap([p for p in fitrange],fit.p)


plt.errorbar(7,c_diff_7ml.mean,xerr=0,yerr=c_diff_7ml.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(5,c_diff_5ml.mean,xerr=0,yerr=c_diff_5ml.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(3,c_diff_3ml.mean,xerr=0,yerr=c_diff_3ml.sdev,fmt='h',mfc='none',color='r',lw=1)

# physical point
plt.errorbar(1, extrapolated[0].mean,xerr=0,yerr=extrapolated[0].sdev,fmt='h',mfc='none', label='extrapolated physical point', color='b', lw=1)
#plt.errorbar(1.05, vc_diff_phys.mean,xerr=0,yerr=vc_diff_phys.sdev,fmt='h',mfc='none', label='actual physical point', color='black', lw=1)

plt.plot([p for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=14,loc='upper right')

plt.xlabel('$m_q/m_l$')
plt.ylabel(r'$\frac{a_{\mu}^{[qcd+qed]}-a_{\mu}^{[qcd]}}{a_{\mu}^{[qcd]}}$', rotation=0, labelpad=50, fontsize=20)

#plt.ylim(bottom=.98)
#plt.xlim(left=0, right=10)

plt.tight_layout()
plt.savefig('../figures/amu_coarsediff_chiExt.png', dpi=500)
plt.show()

plt.close()

