# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re
import gvar as gv


base = '../data/qqed/vcoarse'

name = 'md_vt_vcoarse.gpl'
filepath = os.path.join(base, name)

f = open(filepath, 'r')
f1 = f.readlines()
g = open(base + '/retag_'+ name, 'w+')

retag = ['VTvcmdq0','VTvcmdq1']

tag = gv.dataset.Dataset(filepath).keys()

print('retagging:\n%s'%filepath)
for (a,b) in zip(tag, retag):
    print("%s ----> %s"%(a,b))

for (i,j) in zip(tag,retag):
  for x in f1:  
     y = re.sub(r'\b'+i, j, x)
     if y!=x:
       g.write(y)
    






