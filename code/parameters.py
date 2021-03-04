# PARAMETERS FOR corffitter

TDATA = range(48)
TFIT  = TDATA[4:20]
TP    = 48
NEXP  = range(1,10)
s_coeff = (1, -1)

OSC     = True
CORRFIT = True
PLOT    = False 
SAVEFIG = False


dir = '../data/qed/gpl/'
#dir = '../data/proc_corrs/l3296f211b630m0074m037m440-coul-v5/'
# l3248f211b580m002426m06730m8447
# l3296f211b630m0074m037m440-coul-v5

name = 'm0.0677_Rhox_0.202_chargeAV_outcorr.gpl'
key  = 'RHO_E_AV'    # for fitting, doesn't affect plotting

filename_in = dir + name

corrtag = name[:-4]

otherkey    = None

#priors??

