import dill as pickle
import bz2
import gvar as gv
import sys
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
import g2tools as g2
import matplotlib.pyplot as plt
import corrfitter as cf
import os
from commonweal import base, vt_fname
import numpy as np

## ensemble
a_str = 'vc'

a = (w0/w0overa[a_str])/hbarc           # in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

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

tcut = int(np.ceil(3.5/(a.mean*hbarc)))
tmin = 2
masses = ['3ml', '5ml', '7ml', 'ms']
tcutrange = range(tmin,tcut)

data_unchargedamud = dict.fromkeys(masses)
data_unchargedamuu = dict.fromkeys(masses)
data_chargedamud   = dict.fromkeys(masses)
data_chargedamuu   = dict.fromkeys(masses)

for mq in masses:

    m=mq+':'
    tags = keys[mq]
    print(tags)

    data_unchargedamud[mq] = []
    data_unchargedamuu[mq] = []
    data_chargedamud[mq]   = []
    data_chargedamuu[mq]   = []
    for t in tcutrange:
        vpol = g2.fourier_vacpol(data[tags[0]], Z=ZV, ainv=1/a, tmin=tmin, tmax=t, periodic=True)
        data_unchargedamud[mq].append(g2.a_mu(vpol,1/3.))
        data_unchargedamuu[mq].append(g2.a_mu(vpol,2/3.))

        vpol = g2.fourier_vacpol(data[tags[1]], Z=ZVqedd, ainv=1/a, tmin=tmin, tmax=t, periodic=True)
        data_chargedamud[mq].append(g2.a_mu(vpol,1/3.))

        vpol = g2.fourier_vacpol(data[tags[2]], Z=ZVqedu, ainv=1/a, tmin=tmin, tmax=t, periodic=True)
        data_chargedamuu[mq].append(g2.a_mu(vpol,2/3.))

## find correlations so we can put them back into the model

uncharged_downcorr = gv.evalcorr(data_unchargedamud)
charged_downcorr = gv.evalcorr(data_chargedamud)

## start by doing tcut dependence of the model on its own

model_unchargedamud = dict.fromkeys(masses)
model_unchargedamuu = dict.fromkeys(masses)
model_chargedamud   = dict.fromkeys(masses)
model_chargedamuu   = dict.fromkeys(masses)

# load lsqfit object
fit = pickle.load(bz2.BZ2File('../fits/fitobj/vc_FIT_bin2.pbz2', 'rb'))


#fmt='../figures/plot-{}.png'

#fit.show_plots(save=fmt)


modeluncharged    = gv.BufferDict()
modelchargeddown  = gv.BufferDict()
modelchargedup    = gv.BufferDict()

model_unchargedamud = dict.fromkeys(masses)
model_unchargedamuu = dict.fromkeys(masses)
model_chargedamud   = dict.fromkeys(masses)
model_chargedamuuu   = dict.fromkeys(masses)

for mq in masses:
    m=mq+':'
    tags = ['no-charge', 'down-charge', 'up-charge']

    modeluncharged[mq] = cf.Corr2(datatag=tags[0],tdata=range(tcut),a=(m+'a:n',m+'ao:n'),b=(m+'a:n',m+'ao:n'),dE=(m+'dE:n',m+'dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargeddown[mq] = cf.Corr2(datatag=tags[1],tdata=range(tcut),a=(m+'a:d',m+'ao:d'),b=(m+'a:d',m+'ao:d'),dE=(m+'dE:d',m+'dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargedup[mq] = cf.Corr2(datatag=tags[2],tdata=range(tcut),a=(m+'a:u',m+'ao:u'),b=(m+'a:u',m+'ao:u'),dE=(m+'dE:u',m+'dEo:u'),s=(1.,-1.)).fitfcn(fit.p)
    G = [ gv.corr(g1,g2) for (g1,g2) in zip(modeluncharged[mq],modelchargedup[mq]) ]
    print(G)
    model_unchargedamud[mq] = []
    model_unchargedamuu[mq] = []
    model_chargedamud[mq]   = []
    model_chargedamuu[mq]   = []
    for t in tcutrange:
        vpol = g2.fourier_vacpol(modeluncharged[mq], Z=ZV, ainv=1/a, tmin=tmin, tmax=t, periodic=False)
        model_unchargedamud[mq].append(g2.a_mu(vpol,1/3.))
        model_unchargedamuu[mq].append(g2.a_mu(vpol,2/3.))

        vpol = g2.fourier_vacpol(modelchargeddown[mq], Z=ZVqedd, ainv=1/a, tmin=tmin, tmax=t, periodic=False)
        model_chargedamud[mq].append(g2.a_mu(vpol,1/3.))

        vpol = g2.fourier_vacpol(modelchargedup[mq], Z=ZVqedu, ainv=1/a, tmin=tmin, tmax=t, periodic=False)
        model_chargedamuu[mq].append(g2.a_mu(vpol,2/3.))

model_down_diff = {}
data_down_diff = {}
model_up_diff = {}
data_up_diff = {}


for m in masses:
    model_down_diff[m] = [ (y2-y1) for (y2,y1) in zip(model_chargedamud[m],model_unchargedamud[m]) ] 
    data_down_diff[m] = [ (y2-y1) for (y2,y1) in zip(data_chargedamud[m],data_unchargedamud[m]) ] 
    model_up_diff[m] = [ (y2-y1) for (y2,y1) in zip(model_chargedamuu[m],model_unchargedamuu[m]) ] 
    data_up_diff[m] = [ (y2-y1) for (y2,y1) in zip(data_chargedamuu[m],data_unchargedamuu[m]) ] 


## PLOTTING

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

#fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, figsize=(6,3))

t_fm = [ t*a*hbarc for t in tcutrange ] 

for m in masses[:-3]:
    model_y = model_up_diff[m]
    data_y =  data_up_diff[m]
    plt.errorbar( [t.mean for t in t_fm], [y.mean for y in model_y], xerr=[t.sdev for t in t_fm], yerr=[y.sdev for y in model_y], fmt='h', mfc='none', label=m+' model down diff' )
    plt.errorbar( [t.mean for t in t_fm], [y.mean for y in data_y], xerr=[t.sdev for t in t_fm], yerr=[y.sdev for y in data_y], fmt='h', mfc='none', label=m+' data down diff' )


plt.xlabel('tcut [fm]')
plt.ylabel(r'$\delta a_{\mu}$', rotation=0)
plt.legend()
#plt.savefig('../figures/tcut_dep.png', dpi=500, bbox_inches="tight")
plt.show()


