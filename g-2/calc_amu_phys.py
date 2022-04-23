
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

a_str = 'vcphys' 
B = 1
print("stream [%s] with binsize = %d"%(a_str,B))
f = os.path.join(base, vt_fname[a_str])
a_str = 'vc'
corr = gv.dataset.Dataset(f)

dset = gv.BufferDict()
for tag in corr.keys():
    dset[tag] = 3*np.array(corr[tag])
    print("%d cfgs for tag %s"%(dset[tag].shape[0], tag)) 

corrbin = gv.dataset.bin_data(dset, binsize=B)

data = gv.dataset.avg_data(corrbin)

a = (w0/w0overa[a_str])/hbarc		# in units of (GeV)^-1
ainv = 1/a
print("lattice spacing = %sfm"%(a*hbarc))

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

res = {'diff':[]}

## vc phys analysis

tags = list(data.keys())
print(tags)
vpol = g2.fourier_vacpol(data[tags[0]], Z=ZV, ainv=1/a, periodic=True)
unchargedamuu = g2.a_mu(vpol,2/3.)

vpol = g2.fourier_vacpol(data[tags[1]], Z=ZVqedu, ainv=1/a, periodic=True)
chargedamuu = g2.a_mu(vpol,2/3.)

vpol = g2.fourier_vacpol(data[tags[2]], Z=ZV, ainv=1/a, periodic=True)
unchargedamud = g2.a_mu(vpol,1/3.)

vpol = g2.fourier_vacpol(data[tags[3]], Z=ZVqedd, ainv=1/a, periodic=True)
chargedamud = g2.a_mu(vpol,1/3.)

diff = (chargedamud + chargedamuu) - (unchargedamud + unchargedamuu)
amu  = (unchargedamud + unchargedamuu)
print("results at the physical point on the very coarse ensemble")
print(diff)
print(amu)

inputs = dict(data=data, ZV=ZV, a=a)
outputs = dict(diff=diff, amu=amu)
print(gv.fmt_errorbudget(inputs=inputs, outputs=outputs))



