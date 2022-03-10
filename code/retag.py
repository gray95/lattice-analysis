# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re
import gvar as gv


base = '../data/qqed/ms-tuning-fine/'

name = 'ms0.0368_ps_fine.gpl'
#nname = name[25:]
filepath = os.path.join(base, name)

f = open(filepath, 'r')
f1 = f.readlines()
g = open(base + '/retag_'+ name, 'w+')

retag = ['PSm0368q0','PSm0368q101', 'PSm0368q202']

tag = gv.dataset.Dataset(filepath).keys()

print('retagging:\n%s'%filepath)
for (a,b) in zip(tag, retag):
    print("%s ----> %s"%(a,b))

for i in range(len(tag)):
  for x in f1:  
     y = re.sub(r'\b'+tag[i], retag[i], x)
     if y!=x:
       g.write(y)
    






