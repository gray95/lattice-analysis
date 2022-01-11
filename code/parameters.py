import gvar as gv
import corrfitter as cf
import os
import re
from commonweal import a, hbarc

## flags for both single and matrix fits
OSC       = True
CORRFIT   = True 
PLOT      = False
NOISE     = False 
SAVEFIG   = False 
WRITE_LOG = True 

## MATRIX FIT PARAMETERS
NEXP = range(4,11)            # number of exponentials in fit
t0 = 3                     # initial timeslice to generate priors
tmin = 2
tmax = 6
c_hack = 1									# sometimes -1 needed to generate priors
bin_size = 1

##------------------------------------------##
#l3264f211b600m00507m0507m628a-coul-v5
#l3248f211b580m002426m06730m8447
#l3296f211b630m0074m037m440-coul-v5
#l4864f211b600m001907m05252m6382
ensemble = 'l3296f211b630m0074m037m440-coul-v5/onemm'
corr = 'onemm_fine_m450.gpl'
basedir = '../data/hybrid'
corrpath  = os.path.join(basedir, ensemble, corr)
#SRCs = ['l', 'g']
SRCs = ['R', 'r', 'H', 'h']
KEYFMT = 'onemm.{s1}{s2}'
tag    = re.sub('{s1}{s2}', '', KEYFMT)
ttag   = tag[:-1]
otag   = 'o.'

#dset = gv.dataset.Dataset(corrpath)
data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=tag))     
#data = gv.dataset.avg_data(cf.read_dataset(corrpath))     

s_coeff = (1, -1)
key  = tag+SRCs[0]+SRCs[0]    
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

a_vc    = a[0]					# fermi
a_c     = a[1]
a_f     = a[2]
ainv_vc = (hbarc/a_vc) 			# GeV
ainv_c  = (hbarc/a_c)
ainv_f  = (hbarc/a_f)

