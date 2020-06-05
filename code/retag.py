# A short script to relabel the glasgow correlator data so that the tags are more easily fed 
# into LePage's corrfitter code

import numpy as np
import os
import sys

file = os.path.join('glesga','fine_physical','etacjpsiall.gpl')
corr = np.genfromtxt(file)

np.savetxt('text.txt', corr)

f = open('text.txt', 'r')
f1 = f.readlines()
g = open('etacjpsi-new-new-tags.txt', 'w+')

S = ['l.l', 'l.s1', 'l.s2', 's1.l', 's1.s1', 's1.s2', 's2.l', 's2.s1', 's2.s2']
ETATAGS = ['etac:'+i for i in S]
JPSITAGS = ['jpsi:'+i for i in S]

TAGS = ETATAGS + JPSITAGS

count = 0

for x in f1:  
   x = x.replace('nan', TAGS[count//565])
   g.write(x)
   count=count+1
    

# A more elegant way would be to use re.sub, ie a REGEX





