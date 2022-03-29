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
WRITE_LOG = False 

NEXP = range(1,7)            # number of exponentials in fit
tmin = 3
tmax = 40
bin_size = 1    # implement this

##------------------------------------------##
name = 'ps_fine'
ensemble = 'fine'

corr = name+'.gpl'
basedir = '../../data/qqed'
corrpath  = os.path.join(basedir, ensemble, corr)
tags={'Q0':['PSm0036q0','PSm006q0','PSm0084q0','PSm0364q0'], 'Q101':['PSm0036q101','PSm006q101','PSm0084q101','PSm0364q101'], 'Q202':['PSm0036q202','PSm006q202','PSm0084q202','PSm0364q202']}
#tags={'Q0':['PSm007278q0','PSm01213q0','PSm01698q0','PSm0677q0'], 'Q101':['PSm007278q101','PSm01213q101','PSm01698q101','PSm0677q101'], 'Q202':['PSm007278q202','PSm01213q202','PSm01698q202','PSm0677q202']}
#tags={'Q0':['PSm0362q0','PSm0364q0','PSm0366q0','PSm0368q0'], 'Q101':['PSm0362q101','PSm0364q101','PSm0366q101','PSm0368q101'], 'Q202':['PSm0362q202','PSm0364q202','PSm0366q202','PSm0368q202']}
#tags={ 'Q0':['m00552q0','m0092q0','m01288q0','m0527q0'],'Q101':['m00552q101','m0092q101','m01288q101','m0527q101'],'Q202':['m00552q202','m0092q202','m01288q202','m0527q202'] }
savefile = name+'.p'
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

# LOG
log_name =  'LOG_' + corr        
l = None
if WRITE_LOG:
    l = open(os.path.join('../log', ensemble, log_name), 'w')

a_vc    = w0/w0overa['vc']		# fermi
a_c     = w0/w0overa['c']
a_f     = w0/w0overa['f']
ainv_vc = (hbarc/a_vc) 			# GeV
ainv_c  = (hbarc/a_c)
ainv_f  = (hbarc/a_f)

