import gvar as gv
import corrfitter as cf
import os

##------------------------------------------##
#l3248f211b580m002426m06730m8447\a\
#l3296f211b630m0074m037m440-coul-v5

ensemble = 'l3296f211b630m0074m037m440-coul-v5'
corr = 't0_onemmHy_vec_m0.450.gpl'
basedir = '../data/hybrid'
corrpath  = os.path.join(basedir, ensemble, corr)
SRCs = ['H', 'h']
KEYFMT = 'onemm.{s1}{s2}'
tag    = KEYFMT[:6]
ttag   = tag[:-1]
otag   = ttag + '_o.'

# get temporal extent of lattice
dset = gv.dataset.Dataset(corrpath)
data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=ttag))     
T = data[dset.keys()[0]].size

## single channel parameters
TDATA = range(T)
TFIT  = TDATA[1:20]
TP    = T
NEXP  = range(1,10)
s_coeff = (1,-1)
key  = tag+SRCs[1]+SRCs[1]    
corrtag = corr[:-4]
otherkey    = None

## flags for both single and matrix fits
OSC     = True
CORRFIT = True 
PLOT    = False 
NOISE   = False
SAVEFIG = False
WRITE_LOG = False

# LOG
log_name =  'LOG_' + corr        
l = None
if WRITE_LOG:
    l = open(os.path.join('../log', ensemble, log_name), 'w')


hbarc = 0.197326968
a_vc     = 0.15					# fermi
a_c      = 0.12
a_fine   = 0.09
ainv_vc = (hbarc/a_vc) 			# GeV
ainv_c  = (hbarc/a_c)
ainv_fine=(hbarc/a_fine)

