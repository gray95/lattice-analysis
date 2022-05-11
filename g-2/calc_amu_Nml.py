
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
from commonweal import base, vt_fname

PERIODIC=True
tmax=None
#fit_obj = pickle.load(bz2.BZ2File('./fits/fitobj/vc_FIT_bin2.pbz2', 'rb'))

a_str = 'f' 
B = 2
print("stream [%s] with binsize = %d"%(a_str,B))
f = os.path.join(base, vt_fname[a_str])
corr = gv.dataset.Dataset(f)

dset = gv.BufferDict()
for tag in corr.keys():
    dset[tag] = 3*np.array(corr[tag])
    print("%d cfgs for tag %s"%(dset[tag].shape[0], tag)) 

corrbin = gv.dataset.bin_data(dset, binsize=B)

data = gv.dataset.avg_data(corrbin)

keys = {'3ml':[k for k in data.keys() if '3ml' in k], 
        '5ml':[k for k in data.keys() if '5ml' in k],
        '7ml':[k for k in data.keys() if '7ml' in k],
         'ms':[k for k in data.keys() if 'ms'  in k]}

a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

a_fm = a*hbarc

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

up_noQ = []
up_Q = []
down_noQ = []
down_Q = []

for m in keys:
    tags = keys[m]
    print(tags)

    vpol = g2.fourier_vacpol(data[tags[0]], Z=ZV, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    unchargedamud = g2.a_mu(vpol,1/3.)
    unchargedamuu = g2.a_mu(vpol,2/3.)

    vpol = g2.fourier_vacpol(data[tags[1]], Z=ZVqedd, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    chargedamud = g2.a_mu(vpol,1/3.)

    vpol = g2.fourier_vacpol(data[tags[2]], Z=ZVqedu, tmax=tmax, ainv=1/a, periodic=PERIODIC)
    chargedamuu = g2.a_mu(vpol,2/3.)
    
    amu = chargedamud + chargedamuu
    u_diff = chargedamuu - unchargedamuu
    d_diff = chargedamud - unchargedamud
    diff = u_diff + d_diff

    up_noQ.append(unchargedamuu)
    up_Q.append(chargedamuu)
    down_noQ.append(unchargedamud)
    down_Q.append(chargedamud)

    inputs = dict(data=data, ZV=ZVqedu, a=a)
    outputs = dict(updiff=u_diff)
    print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs, verify=True))

res = gv.BufferDict()
res['up-nocharge'] = up_noQ
res['down-nocharge'] = down_noQ
res['up-charge'] = up_Q
res['down-charge'] = down_Q

## save results to disk
#gv.dump( res, 'out/'+a_str+'_no_vector_fit.p', add_dependencies=True)

