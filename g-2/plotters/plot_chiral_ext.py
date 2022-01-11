import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
from commonweal import vc_rt, c_rt , c_diff15, c_diff, vc_diff, vc_diff15
from gvar import log
import sys

ydata = vc_diff15
xdata = np.array([3,5,7])

def extrap(x,p):
    return ( p[0] + p[1]*x )   
def fitargs(z):
    dp0 = z
    dp1 = z
    prior = gv.gvar([gv.gvar(-1.5e-11, dp1), gv.gvar(0, dp1)])
    return dict(prior=prior, fcn=extrap, data=(xdata, ydata))

#ef extrap(x,p):
#   res = []
#   for i,point in enumerate(x):
#       res.append(p['a'][0] * (1 + p['a'][1]*(point) )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
#   return res

#def fitargs(z):
#    dp = z  #z['dp0']
#    #dpp = z['dp1']
#    prior = {}
#    prior['a'] = gv.gvar([gv.gvar(0.0, dp), gv.gvar(0.2, dp)])
#    return dict(prior=prior, fcn=extrap, data=([7, 5, 3], list(reversed(c_diff))))


hbarc = 0.197326968

## EMPIRICAL BAYES
#z0 = { 'dp0' : 0.1, 'dp1' : '0.0' }
z0 = 0.1
fit, z = lsq.empbayes_fit(z0, fitargs)
print(fit.format(True))
print("prior width that maxmises logGBF: ", z)
#print("prior width on a[1] that maxmises logGBF: ", z['dp1'])

#prior = {}
#prior['a'] = [gv.gvar('1.0(1)'),gv.gvar('0.0(1)')] #,gv.gvar(0,1),gv.gvar(0,1)]
## extrapolation
#fit = lsq.nonlinear_fit(data=([3,5,7], ydata), prior=prior, fcn=extrap)
#print("Fit results")
#print(fit)

#extrapolated = extrap([1], fit.p)
extrapolated = extrap(1, fit.p)
print("physical value: ", extrapolated)

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
#plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')


fitrange = np.arange(0,8,0.01)
fitline = [extrap(x,fit.p) for x in fitrange]

plt.errorbar(7,ydata[2].mean,xerr=0,yerr=ydata[2].sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(5,ydata[1].mean,xerr=0,yerr=ydata[1].sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(3,ydata[0].mean,xerr=0,yerr=ydata[0].sdev,fmt='h',mfc='none',color='r',lw=1)

# physical point
plt.errorbar(1.05, extrapolated.mean,xerr=0,yerr=extrapolated.sdev,fmt='h',mfc='none', label='linear extrapolation', color='b', lw=1)


plt.plot([p for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=14,loc='upper right')

plt.xlabel(r'$\frac{m_q}{m_l}$', fontsize=25, labelpad=20)
#plt.ylabel(r'$\frac{a_{\mu}^{\mathrm{qcd+qed}}}{a_{\mu}^{\mathrm{qcd}}}$', rotation=0, labelpad=25, fontsize=20)
plt.ylabel(r'$\delta a_{\mu}$', rotation=0, labelpad=25, fontsize=20)
plt.title('0.12fm ensemble - 694 configs')

#plt.ylim(top=0.00)
#plt.xlim(left=0, right=10)

plt.tight_layout()
#plt.savefig('../figures/amu_c_diff_15chiXt.png', dpi=500)
plt.show()

plt.close()

