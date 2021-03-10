import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt

a_vcoarse = (gv.gvar('0.1715(9)')/gv.gvar('1.13215(35)'))
a_coarse = (gv.gvar('0.1715(9)')/gv.gvar('1.41490(60)'))
a_fine = (gv.gvar('0.1715(9)')/gv.gvar('1.9518(7)'))

# data hard wired
vc = gv.gvar('3.30(43)')
c = gv.gvar('4.60(59)')
f = gv.gvar('8.4(1.8)')

# finite volume  Date: 2 June 2020 at 21:51:12 BST, email from Peter
vcfv = gv.gvar('6.8(1.4)')
cfv = gv.gvar('5.4(1.1)')
ffv = gv.gvar('3.58(72)')

##  correlation
vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,0.5,0.5],[0.5,1,0.5],[0.5,0.5,1]])
#vcfv,cfv,ffv = gv.correlate([vcfv,cfv,ffv],[[1,1,1],[1,1,1],[1,1,1]])

vccorr = vc+vcfv
ccorr = c+cfv
fcorr = f+ffv

print(gv.evalcorr([vccorr,ccorr,fcorr]))

def extrap(x,p):
    res = []
    for i,point in enumerate(x):
        res.append(p['a'][0] * (1 + p['a'][1]*(point)**2 )) # + p['a'][2]*(point*0.5)**4 )) # + p['a'][3]*mistunings[i]))
    return res

hbarc = 0.197326968

prior = {}
prior['a'] = [gv.gvar('12(10)'),gv.gvar(0,100),gv.gvar(0,1),gv.gvar(0,1)]

print(1/(a_vcoarse.mean/hbarc),1/(a_coarse.mean/hbarc))

## a**2 extrapolation
fit = lsq.nonlinear_fit(data=([a_vcoarse.mean/hbarc,a_coarse.mean/hbarc,a_fine.mean/hbarc],[vccorr,ccorr,fcorr]),prior=prior,fcn=extrap)
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

plt.errorbar(a_vcoarse.mean**2,vc.mean,xerr=0,yerr=vc.sdev,fmt='h',mfc='none',color='r',label='Fermilab/HPQCD/MILC',lw=1)
plt.errorbar(a_coarse.mean**2,c.mean,xerr=0,yerr=c.sdev,fmt='h',mfc='none',color='r',lw=1)
plt.errorbar(a_fine.mean**2,f.mean,xerr=0,yerr=f.sdev,fmt='h',mfc='none',color='r',lw=1)

plt.errorbar(a_vcoarse.mean**2,vccorr.mean,xerr=0,yerr=vccorr.sdev,fmt='h',color='r',label='Fermilab/HPQCD/MILC corrected',lw=1)
plt.errorbar(a_coarse.mean**2,ccorr.mean,xerr=0,yerr=ccorr.sdev,fmt='h',color='r',lw=1)
plt.errorbar(a_fine.mean**2+0.0005,fcorr.mean,xerr=0,yerr=fcorr.sdev,fmt='h',color='r',lw=1)

plt.plot([p**2 for p in fitrange],[p.mean for p in fitline],'--',color='r')
plt.fill_between([p**2 for p in fitrange],[p.mean-p.sdev for p in fitline],[p.mean+p.sdev for p in fitline],alpha=0.5,lw=0,color='r')

arbc = 0.11406

plt.errorbar(arbc**2,11.2,xerr=0,yerr=4,fmt='s',mfc='none',color='g',label='RBC/UKQCD [1801.07224]',lw=1)

#plt.errorbar(0,12.8,xerr=0,yerr=1.9,fmt='s',color='b',label='BMW [1711.04980]',lw=1)
plt.errorbar(0,13.15,xerr=0,yerr=1.8,fmt='s',color='b',label='BMW [2002.12347]',lw=1)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
plt.legend(handles=handles,labels=labels,frameon=False,fontsize=8,loc='lower left')

plt.xlabel('$a^2\ [\mathrm{fm}^2]$')
plt.ylabel('$-a_{\mu}^{\mathrm{disconn.}} \\times 10^{10}$')

plt.ylim(bottom=-5)
plt.xlim(left=-0.0005)

plt.tight_layout()
#plt.savefig('./amu_disconn_summary.png')
plt.show()

plt.close()

