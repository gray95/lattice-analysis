# PARAMETERS FOR corffitter

T     = 48
TDATA = range(T)
TFIT  = TDATA[8:16]
TP    = T
NEXP  = range(1,10)
s_coeff = (-1, 1)

OSC     = True
CORRFIT = True
PLOT    = False 
SAVEFIG = False


#dir = '../data/qed/gpl_old/'
#dir = '../data/proc_corrs/l3296f211b630m0074m037m440-coul-v5/'
dir  = '../data/proc_corrs/l3248f211b580m002426m06730m8447/a/'
# l3296f211b630m0074m037m440-coul-v5

name = 't0_onemp_ranwall_m0.8447.txt'
key  = 'onemp.gl'    # for fitting, doesn't affect plotting

filename_in = dir + name

corrtag = name[:-4]

otherkey    = None


