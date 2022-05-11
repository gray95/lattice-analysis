import gvar as gv
import corrfitter as cf
import os
import re
import sys
sys.path.append('../plotters')
from commonweal import w0overa, hbarc, w0

CORRFIT   = True 
PLOT      = False
NOISE     = False 
SAVEFIG   = False 

LBL = 'vc'
NEXP = range(1,7)            # number of exponentials in fit
tmin = 3
tmax = 20
bnsze = 1     
norm = 3
##------------------------------------------##

if LBL=='vc':
    name = 'ps_vcoarse'
    ensemble = 'vcoarse'
elif LBL=='c':
    name = 'ps_coarse'
    ensemble = 'coarse'
elif LBL=='f':
    name = 'ps_fine'
    ensemble = 'fine'

corr = name+'.gpl'
basedir = '../../data/qqed'
corrpath  = os.path.join(basedir, ensemble, corr)
TAGS={'f':{'Q0':['PSf3mlq0','PSf5mlq0','PSf7mlq0','PSfmsq0'], 'Q101':['PSf3mlq1','PSf5mlq1','PSf7mlq1','PSfmsq1'], 'Q202':['PSf3mlq2','PSf5mlq2','PSf7mlq2','PSfmsq2']}, 'c':{'Q0':['PSc3mlq0','PSc5mlq0','PSc7mlq0','PScmsq0'], 'Q101':['PSc3mlq1','PSc5mlq1','PSc7mlq1','PScmsq1'], 'Q202':['PSc3mlq2','PSc5mlq2','PSc7mlq2','PScmsq2']}, 'vc':{'Q0':['PSvc3mlq0','PSvc5mlq0','PSvc7mlq0','PSvcmsq0'], 'Q101':['PSvc3mlq1','PSvc5mlq1','PSvc7mlq1','PSvcmsq1'], 'Q202':['PSvc3mlq2','PSvc5mlq2','PSvc7mlq2','PSvcmsq2']}}
tags = TAGS[LBL]
savefile = './fits/'+name+'.p'
NMASSES = len(tags['Q0'])
data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=tags['Q0'][0]))     

key  = tags['Q0'][0]    
corrtag = corr[:-4]
T = data[key].size
TDATA = range(T)
TP    = T
TFIT  = TDATA[tmin:tmax+1]

# print no of cfgs
print(corrpath)
print("time extent: %d"%T)
print("fit range: [%d,%d]"%(tmin,tmax))
print("fitting %d quark masses"%NMASSES)

a    = w0/w0overa[LBL]		# fermi
ainv = (hbarc/a) 		# GeV

