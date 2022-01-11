import gvar as gv
import corrfitter as cf
import os
import re
import sys
sys.path.append('plotters')
from commonweal import w0overa, hbarc, w0

OSC       = False
CORRFIT   = True 
PLOT      = True
NOISE     = False 
SAVEFIG   = False 
WRITE_LOG = False 

NEXP = range(1,10)            # number of exponentials in fit
tmin = 2
tmax = 20
bin_size = 1

##------------------------------------------##
ensemble = 'coarse'
corr = 'm0.0527_Pseuodoscalar_Q0.202_chargeAV_outcorr.gpl'
basedir = '../data/qqed'
corrpath  = os.path.join(basedir, ensemble, corr)
#KEYFMT = 'onemm.{s1}{s2}'
#tag    = re.sub('{s1}{s2}', '', KEYFMT)
#ttag   = tag[:-1]
#otag   = 'o.'
tag='PSEUDO_E_AV_Q0.202'

data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=tag))     

s_coeff = (1, -1)
key  = tag    
corrtag = corr[:-4]
#otherkey = tag+SRCs[0]+SRCs[0] 
otherkey = None
T = data[key].size
TDATA = range(T)
TP    = T
TFIT  = TDATA[tmin:tmax+1]

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

