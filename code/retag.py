# A short script to relabel correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import os
import sys
import re

name = 'PION_PS_B_I_c_d_m0.8447b.txt'
filepath = os.path.join('../data/proc_corrs/l3248f211b580m002426m06730m8447', name)


f = open(filepath, 'r')
f1 = f.readlines()
g = open('retag_'+name, 'w+')

for x in f1:  
   x = re.sub('PION_PS_B_I_c_d_m0.8447', 'onemp_0_src_m0.8447b', x)
   g.write(x)
    






