import dill as pickle
import bz2
import gvar as gv
import sys
sys.path.append('..')
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
import g2tools as g2
import matplotlib.pyplot as plt
import corrfitter as cf
import os
from commonweal import base, vt_fname
import numpy as np

## ensemble
a_str = 'vc'
# load lsqfit object
fit= pickle.load(bz2.BZ2File('../fits/fitobj/vtFIT_vc_bin2_tmin2.pbz2', 'rb'))
a = (w0/w0overa[a_str])/hbarc           # in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))
#tcut = int(np.ceil(6./(a.mean*hbarc)))
tmin = 2
masses = ['3ml', '5ml', '7ml', 'ms']


ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

B = 2
print("stream [%s] with binsize = %d"%(a_str,B))
f = os.path.join(base, vt_fname[a_str])
corr = gv.dataset.Dataset(f)

dset = gv.BufferDict()

for tag in corr.keys():
    corr[tag] = 3*np.array(corr[tag])


corrbin = gv.dataset.bin_data(corr, binsize=B)
print("%d cfgs"%corrbin.samplesize)

data = gv.dataset.avg_data(corrbin)

keys = {'3ml':[k for k in data.keys() if '3ml' in k], 
        '5ml':[k for k in data.keys() if '5ml' in k],
        '7ml':[k for k in data.keys() if '7ml' in k],
         'ms':[k for k in data.keys() if 'ms'  in k]}

print(keys)

T = data[keys['3ml'][0]].shape[0]

tleft = 2*T//3
tright = T

data_t_fm = [ t*a*hbarc for t in range(T//2-5,T//2) ] 
model_t_fm = [ t*a*hbarc for t in range(tleft,tright) ] 

modeluncharged    = gv.BufferDict()
modelchargeddown  = gv.BufferDict()
modelchargedup    = gv.BufferDict()

for mq in masses:
    m=mq+':'
    tags = ['no-charge', 'down-charge', 'up-charge']

    modeluncharged[mq] = cf.Corr2(datatag=tags[0],tdata=range(T),a=(m+'a:n',m+'ao:n'),b=(m+'a:n',m+'ao:n'),dE=(m+'dE:n',m+'dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargeddown[mq] = cf.Corr2(datatag=tags[1],tdata=range(T),a=(m+'a:d',m+'ao:d'),b=(m+'a:d',m+'ao:d'),dE=(m+'dE:d',m+'dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargedup[mq] = cf.Corr2(datatag=tags[2],tdata=range(T),a=(m+'a:u',m+'ao:u'),b=(m+'a:u',m+'ao:u'),dE=(m+'dE:u',m+'dEo:u'),s=(1.,-1.)).fitfcn(fit.p)

## PLOTTING

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

#fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, figsize=(6,3))

for m in masses:
    print("plotting mass %s"%m)
    #plt.clf()
    data_y =  data[keys[m][0]][T//2-5:T//2]
    model_y = modeluncharged[m][tleft:tright]
    plt.errorbar( [t.mean for t in model_t_fm], [y.mean for y in model_y], xerr=[t.sdev for t in model_t_fm], yerr=[y.sdev for y in model_y], fmt='h', mfc='none', alpha=0.6, color='red', label=m+' model' )
    #plt.errorbar( [t.mean for t in data_t_fm], [y.mean for y in data_y], xerr=[t.sdev for t in data_t_fm], yerr=[y.sdev for y in data_y], fmt='h', alpha=0.6, color='black', label=m+' data' )
    #plt.errorbar( [t.mean for t in model_t_fm], [abs(y.sdev/y.mean) for y in model_y], fmt='h', mfc='none', alpha=0.6, color='red', label=m+' model frac error' )
    #plt.errorbar( [t.mean for t in data_t_fm], [abs(y.sdev/y.mean) for y in data_y], fmt='h', alpha=0.6, color='black', label=m+' data frac error' )
    plt.xlabel('t [fm]')
    plt.axhline(y=0,linewidth=0.5, linestyle='--')
    #plt.ylabel(r'$G(t)$', rotation=0, labelpad=30, fontsize=16)
    plt.legend()
    plt.show()
    #plt.savefig('../figures/tcutdep_zoomup_'+a_str+'_'+m+'.png', dpi=500, bbox_inches="tight")


