# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re
import gvar as gv
from commonweal import onemp, onemp_paths

def retag_corr(path_to_corr, newtags):

    f = open(path_to_corr, 'r')
    f1 = f.readlines()

    d, corr = os.path.split(path_to_corr) 
    g = open(d + '/retag_'+ corr, 'w+')

    retag = newtags
    tag = gv.dataset.Dataset(path_to_corr).keys()

    print('retagging:\n%s'%path_to_corr)
    for (a,b) in zip(tag, retag):
        print("%s ----> %s"%(a,b))

    for (i,j) in zip(tag,retag):
      for x in f1:  
#          y = re.split('\s+', x)[0]
#          if y==i:
#              g.write(x)
         y = re.sub(r'\b'+i, j, x)
         if y!=x:
             g.write(y)
     
## g2
#base = '../data/qqed/fine'
#masses = [0.0036, 0.006, 0.0084, 0.0364]
#m = ['3ml','5ml','7ml','ms']
#for M in zip(masses,m):
#    ps_filepath = 'm'+str(M[0])+'_ps_fine.gpl'
#    ps_retag = ['PSf'+M[1]+'q0', 'PSf'+M[1]+'q1','PSf'+M[1]+'q2']    
#    retag_corr(base, ps_filepath, ps_retag) 
#
#    vt_filepath = 'm'+str(M[0])+'_vt_fine.gpl'
#    vt_retag = ['VTf'+M[1]+'q0', 'VTf'+M[1]+'q1','VTf'+M[1]+'q2']    
#    retag_corr(base, vt_filepath, vt_retag) 

##hybrids
newtags = ['onemp.ll', 'onemp.gl', 'onemp.lg', 'onemp.gg']
retag_corr(onemp_paths['c'], newtags)




