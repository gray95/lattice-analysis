# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re

base = '../data/qed/vcoarse/dan_b'

name = 'm0.001524_rhox_0.202cav_outcorr.gpl'
nname = name[25:]
filepath = os.path.join(base, name)

f = open(filepath, 'r')
f1 = f.readlines()
g = open(base + '/retag_'+ name, 'w+')


for x in f1:  
   x = re.sub('0.001524_GXGX', 'rho_m0.001524_cav', x)
   g.write(x)
    






