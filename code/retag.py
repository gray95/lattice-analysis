# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re
import gvar as gv


base = '../data/qqed/coarse/'

name = '3ml_rho_coarse.gpl'
#nname = name[25:]
filepath = os.path.join(base, name)

f = open(filepath, 'r')
f1 = f.readlines()
g = open(base + '/retag_'+ name, 'w+')

retag = ['VTc3mlq0','VTc3mlq2', 'VTc3mlq1']

tag = gv.dataset.Dataset(filepath).keys()

print('retagging:\n%s'%filepath)
for (a,b) in zip(tag, retag):
    print("%s ----> %s"%(a,b))

for (i,j) in zip(tag,retag):
  for x in f1:  
     y = re.sub(r'\b'+i, j, x)
     if y!=x:
       g.write(y)
    






