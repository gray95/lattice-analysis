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

# data hard wired
svn = gv.gvar('0.990(23)')
fve = gv.gvar('0.988(40)')
thr = gv.gvar('0.969(73)')
rphys = gv.gvar('1.017(84)')

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('12(10)'),gv.gvar(0,100),gv.gvar(0,1),gv.gvar(0,1)]


## extrapolation
fit = lsq.nonlinear_fit(data=([seven_ml/hbarc,five_ml/hbarc,three_ml/hbarc],[svn,fve,thr]),prior=prior,fcn=extrap)
print("Fit results")
print(fit)

print("physical value of anomaly ", extrap([mphys], fit.p))

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')

fitrange = np.arange(0,seven_ml,0.0001)
fitline = extrap([p/hbarc for p in fitrange],fit.p)

plt.errorbar(seven_ml,svn.mean,xerr=0,yerr=svn.sdev,fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)
plt.errorbar(five_ml,fve.mean,xerr=0,yerr=fve.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(three_ml,thr.mean,xerr=0,yerr=thr.sdev,fmt='h',mfc='none',color='r',lw=1)

# physical point
plt.errorbar(mphys/mphys,rphys.mean,xerr=0,yerr=rphys.sdev,fmt='h',mfc='none', label='physical point', color='b', lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower left')

plt.xlabel('$(m_q/m_l)$')
plt.ylabel(r'$\frac{a_{\mu}^{\mbox{QCD+QED}}}{a_{\mu}^{\mbox{QCD}}}$')

plt.ylim(bottom=.8)
plt.xlim(left=0, right=10)

plt.tight_layout()
plt.savefig('amu_ratio_chiExt.png', dpi=600)
plt.show()

plt.close()

