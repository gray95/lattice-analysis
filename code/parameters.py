# PARAMETERS FOR corffitter

TDATA = range(48)
TFIT  = TDATA[2:16]
TP    = 48
NEXP  = range(1,7)
s_coeff = (-1, -1)

OSC     = True
CORRFIT = True
PLOT    = True
SAVEFIG = False

dir = 'proc_corrs/ape4/'
# l3248f211b580m002426m06730m8447

name = 'hybrid_lrho_4ape.txt'
key         = '1-+_4ape_lrho'    # for fitting, doesn't affect plotting

filename_in = dir + name

corrtag = name[:-4]

otherkey    = None

#priors??

