import gvar as gv
import corrfitter as cf
import os
import re
from commonweal import w0, w0overa, hbarc
from commonweal import ensembles, onemp_paths, onemm_paths

## flags for both single and matrix fits
OSC       = True
CORRFIT   = True 
NOISE     = True 
EIG       = False

PLOT      = False
SAVEFIG   = False

SAVE_RES  = False
WRITE_LOG = False

# which ensemble [ vc c c-10 f f-5 ]
e_str = 'f-5'    

## MATRIX FIT PARAMETERS
NEXP = range(1,10)            # number of exponentials in fit
t0 = 3                    # initial timeslice to generate priors
tmin = 2
tmax = 6 
D = 6                   # off-diagonal fitting range = TFIT[:D]
c_hack = 1		# sometimes -1 needed to generate priors
bin_size = 2 
s_coeff = (1,-1)

corrpath  = onemm_paths[e_str]
#SRCs      = ['l','g']
SRCs = ['H', 'h', 'R', 'r']
EIG_SRCs  = ['0','1', '2', '3']
KEYFMT = 'onemm.{s1}{s2}'
key  = re.sub('{s1}{s2}', SRCs[0]+SRCs[0], KEYFMT)    
otherkey = None
tag    = '1mm.'
otag   = 'pp.'

## save masses and amplitudes in gvar object
bsedir = '/home/gray/Desktop/lattice-analysis/code/out/'
save_res_name = bsedir+e_str+'/onemp.p'
save_fit      = 'out/'+e_str+'/onemp.pbz2'

# LOG
log_name =  'LOG_' + tag + e_str        
l = None
if WRITE_LOG:
    l = open(os.path.join('../log', ensembles[e_str], log_name), 'w')
    print("writing to logfile %s"%l)
