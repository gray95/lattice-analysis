#!/usr/local/bin/python3


import math

corr_tag = "RHO_E_AV"

##ainv_gev = 0.1973 / 0.15 
ainv_gev = 1.0

#  now set in  plot_corr.py 
nrm = 1
##nrm = 3.0*32*32*32

##nrm = 32*32*32*9

##
##
##

import sys
import numpy as np
import math
import matplotlib.pyplot as plt


# my modules
import sys
sys.path.append('/Users/cmcneile/projects/hybrids/fine/1-+/my_fit/src')

import util

##mass = "m0.001524"
mass = "m0.0677"
#  filename with the correlators
##fff = "Pseuodoscalar_0.202_chargeAV_outcorr.gpl"
fff = "gpl/" + mass  +  "_Rhox_0.202_chargeAV_outcorr.gpl"

nt = 48
no_config = 101


##
##  new data
##

corr = util.load_data(fff, nt,no_config, corr_tag  )
##corr *= -1 
corr /= nrm



print("Normalizatio factor " , nrm, " applied")
print("Computing jackknife correlators")
tt, corr_mean, corr_err = util.calc_corr(corr, nt,no_config, 0.0) 

## 
##  Dan's data
##
nconfig_dan = 381

corr_dan  = util.load_data("../docs/rho_vcphys_bothcharges_m0.001524.gpl", nt, nconfig_dan, "charged-up"  )
tt_dan, corr_dan_mean, corr_dan_err = util.calc_corr(corr_dan, nt,nconfig_dan, 0.1) 

##
##  setup plots
##

##plt.title("Pseuodo-scalar meson correlator")
plt.title("Vector meson correlator (up quark mass )")
##plt.legend()

tt *= 0.15 
tt_dan *= 0.15 


plt.errorbar(tt, corr_mean , corr_err  ,  fmt= 'ro' , label="new meas (stream b)  (nstats = 101)" )
plt.errorbar(tt_dan, corr_dan_mean , corr_dan_err  ,  fmt= 'bo' , label="Dan meas. (stream a) (nstats = 381)" )

plt.yscale('log')
plt.xlabel('t fm')
plt.ylabel('corr(t)')
plt.xlim(0,2)
##plt.xlim(0,0.25)



##plt.ylim(0,8)
plt.legend()

plt.savefig("figures/rho_corr.png")
plt.show()

