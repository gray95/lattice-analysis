import gvar as gv
import corrfitter as cf
import os

##------------------------------------------##
ensemble = 'l3296f211b630m0074m037m440-coul-v5'
corr = 'comb8_onempHy_m0.450.txt'
basedir = '../data/hybrid'
corrpath  = os.path.join(basedir, ensemble, corr)

KEYFMT = 'onemp.{s1}{s2}'
tag    = KEYFMT[:6]
ttag   = tag[:-1]
otag   = ttag + '_o.'

# get temporal extent of lattice
dset = gv.dataset.Dataset(corrpath)
data = gv.dataset.avg_data(cf.read_dataset(corrpath, grep=ttag))     
T = data[dset.keys()[0]].size

## single channel parameters
TDATA = range(T)
TFIT  = TDATA[3:15]
TP    = T
NEXP  = range(1,10)
s_coeff = (1, -1)
key  = 'onemp.ll'    
corrtag = corr[:-4]
otherkey    = None

## flags for both single and matrix fits
OSC     = True
CORRFIT = True
PLOT    = True 
NOISE   = False
SAVEFIG = False
WRITE_LOG = False

# LOG
log_name =  'LOG_' + corr        
l = None
if WRITE_LOG:
    l = open(os.path.join('../log', ensemble, log_name), 'w')


hbarc = 0.197326968
a     = 0.15					# fermi
ainv = (hbarc/a) 			# GeV
