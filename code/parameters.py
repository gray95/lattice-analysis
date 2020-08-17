# PARAMETERS FOR corffitter

TDATA = range(96)
TFIT  = TDATA[4:20]
TP    = 96
NEXP  = range(1,10)
s_coeff = (1, 1)

OSC     = True
CORRFIT = False
PLOT    = True 
SAVEFIG = False

dir = '../data/proc_corrs/l3296f211b630m0074m037m440-coul-v5-new/'
# l3248f211b580m002426m06730m8447
# l3296f211b630m0074m037m440-coul-v5

name = 't0_onempHy_m0.450.txt'
key  = 'onemp.ll'    # for fitting, doesn't affect plotting

filename_in = dir + name

corrtag = name[:-4]

otherkey    = None

#priors??

