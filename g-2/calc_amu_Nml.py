
import os
import sys
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

a_str = 'vc' 
B = 1
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

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

res = gv.BufferDict()
res['up-nocharge'] = {}
res['up-charge']   = {}
res['up-diff']     = {}
res['down-nocharge'] = {}
res['down-charge']   = {}
res['down-diff']     = {}

for m in keys:
    tags = keys[m]
    print(tags)
    vpol = g2.fourier_vacpol(data[tags[0]], Z=ZV, ainv=1/a, periodic=True)
    unchargedamud = g2.a_mu(vpol,1/3.)
    unchargedamuu = g2.a_mu(vpol,2/3.)

    vpol = g2.fourier_vacpol(data[tags[1]], Z=ZVqedd, ainv=1/a, periodic=True)
    chargedamud = g2.a_mu(vpol,1/3.)

    vpol = g2.fourier_vacpol(data[tags[2]], Z=ZVqedu, ainv=1/a, periodic=True)
    chargedamuu = g2.a_mu(vpol,2/3.)
    
    if m=='ms':
        diff = chargedamud - unchargedamud
        amu = chargedamud
        res['down-diff'][m] = diff
        res['down-charge'][m] = chargedamud
        res['down-nocharge'][m] = unchargedamud
    else:
        amu = chargedamud + chargedamuu
        u_diff = chargedamuu - unchargedamuu
        d_diff = chargedamud - unchargedamud
        diff = u_diff + d_diff
        res['up-diff'][m] = u_diff
        res['down-diff'][m] = d_diff
        res['up-charge'][m] = chargedamuu
        res['up-nocharge'][m] = unchargedamuu
        res['down-charge'][m] = chargedamud
        res['down-nocharge'][m] = unchargedamud

    print("################ %s results ################"%m)
    print(diff)
    print(amu)

    inputs = dict(data=data, ZV=ZV, a=a)
    outputs = dict(diff=diff)
    print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs))

## save results to disk
#gv.dump( res, 'out/'+a_str+'_no_vector_fit.p', add_dependencies=True)
