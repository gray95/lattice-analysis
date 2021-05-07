# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re
import gvar as gv


base = '../data/qed/vcoarse/dan_a'

name = 'rho_vcphys_bothcharges_up_down.gpl'
#nname = name[25:]
filepath = os.path.join(base, name)

f = open(filepath, 'r')
f1 = f.readlines()
g = open(base + '/retag_'+ name, 'w+')

retag = ['rho_m0.001524', 'rho_m0.001524_ucav', 'rho_m0.003328', 'rho_m0.003328_dcav']

tag = gv.dataset.Dataset(filepath).keys()
print(tag, ' getting replaced with: ', retag)

for i in range(4):
  for x in f1:  
     y = re.sub(r'\b'+tag[i], retag[i], x)
     if y!=x:
       g.write(y)
    






