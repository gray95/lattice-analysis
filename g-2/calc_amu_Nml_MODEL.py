
import os
import sys
import dill as pickle
import bz2
sys.path.append('./plotters')
import corrfitter as cf
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import g2tools as g2
import lsqfit as lsq
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
from fitting import make_data
from commonweal import base, vt_fname, vtmodels, vtfitobj
import blinding
from corrfitter import Corr2, CorrFitter


PERIODIC=False
tmax=None
tcut = 100

##blinding
bF = "g2-conn-blindfactor.bin"
print("load blinding factor from %s"%bF)
blindF = blinding.read_blind(bF)

a_str = 'f' 

B = 2
print("stream [%s] with binsize = %d"%(a_str,B))
f = os.path.join(base, vt_fname[a_str])
corr = gv.dataset.Dataset(f)

dset = gv.BufferDict()

for tag in corr.keys():
    corr[tag] = 3*blindF*np.array(corr[tag])


corrbin = gv.dataset.bin_data(corr, binsize=B)
print("%d cfgs"%corrbin.samplesize)

data = gv.dataset.avg_data(corrbin)

keys = {'3ml':[k for k in data.keys() if '3ml' in k], 
        '5ml':[k for k in data.keys() if '5ml' in k],
        '7ml':[k for k in data.keys() if '7ml' in k],
         'ms':[k for k in data.keys() if 'ms'  in k]}

print(keys)

fit = pickle.load(bz2.BZ2File(vtfitobj[a_str], 'rb'))

a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

a_fm = a*hbarc
tstar_fm = 1.5
tstar = int(np.ceil(tstar_fm/(a.mean*hbarc)))
print("replace data with model after %f fm (t=%d)"%(tstar_fm,tstar))

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

up_noQ = []
up_Q = []
down_noQ = []
down_Q = []

newdata=gv.BufferDict()

for m in keys:
    tags=keys[m]
    print(tags)    
    fitdatauncharged = blindF*Corr2(datatag=tags[0],tdata=range(tcut),a=(m+':a:n',m+':ao:n'),b=(m+':a:n',m+':ao:n'),dE=(m+':dE:n',m+':dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
    fitdatachargeddown = blindF*Corr2(datatag=tags[1],tdata=range(tcut),a=(m+':a:d',m+':ao:d'),b=(m+':a:d',m+':ao:d'),dE=(m+':dE:d',m+':dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
    fitdatachargedup = blindF*Corr2(datatag=tags[2],tdata=range(tcut),a=(m+':a:u',m+':ao:u'),b=(m+':a:u',m+':ao:u'),dE=(m+':dE:u',m+':dEo:u'),s=(1.,-1.)).fitfcn(fit.p)
    
    newdata[tags[0]] = np.append( data[tags[0]][:tstar+1], fitdatauncharged[tstar+1:] )
    newdata[tags[1]] = np.append( data[tags[1]][:tstar+1], fitdatachargeddown[tstar+1:] )
    newdata[tags[2]] = np.append( data[tags[2]][:tstar+1], fitdatachargedup[tstar+1:] )

    vpol = g2.fourier_vacpol(newdata[tags[0]], Z=ZV, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    unchargedamud = g2.a_mu(vpol,1/3.)
    unchargedamuu = g2.a_mu(vpol,2/3.)

    vpol = g2.fourier_vacpol(newdata[tags[1]], Z=ZVqedd, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    chargedamud = g2.a_mu(vpol,1/3.)

    vpol = g2.fourier_vacpol(newdata[tags[2]], Z=ZVqedu, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    chargedamuu = g2.a_mu(vpol,2/3.)
    
    amu = chargedamud + chargedamuu
    u_diff = chargedamuu - unchargedamuu
    d_diff = chargedamud - unchargedamud
    diff = u_diff + d_diff
    if m=='ms':
        amu = chargedamud
        diff = chargedamud - unchargedamud
 
    print(amu)
    print(diff)
    up_noQ.append(unchargedamuu)
    up_Q.append(chargedamuu)
    down_noQ.append(unchargedamud)
    down_Q.append(chargedamud)

    inputs = dict(data=data, ZV=ZVqedu, a=a)
    outputs = dict(updiff=u_diff)
    #print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs, verify=True))

res = gv.BufferDict()
res['up-nocharge'] = up_noQ
res['down-nocharge'] = down_noQ
res['up-charge'] = up_Q
res['down-charge'] = down_Q

## save results to disk
#gv.dump( res, 'out/'+a_str+'_nofit.p', add_dependencies=True)

