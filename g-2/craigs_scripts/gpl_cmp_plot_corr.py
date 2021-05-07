#!/Users/cmcneile/anaconda3/bin/python3.6

#
#  compare the positive and negative
#

import math
import sys

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
nt = 48

#  filename with the correlators
##fff = "Pseuodoscalar_0.202_chargeAV_outcorr.gpl"


def plot_cxorr(fff_, no_config_, fmt_, lab_, sh, corr_tag) :
  corr = util.load_data(fff, nt,no_config_, corr_tag)

  print("Computing jackknife correlators")
  tt, corr_mean, corr_err = util.calc_corr(corr, nt,no_config_,sh ) 
  plt.errorbar(tt, corr_mean , corr_err  ,  fmt= fmt_ , label=lab_ )


#
#  
#

if 1 :
  fff = "gpl/m0.001524_Rhox_0.202_chargeAV_outcorr.gpl"
  no_config = 184
  plot_cxorr(fff, no_config, "ro","Charge average  (nstats = 101)" , 0.0, "RHO_E_AV") 

if 1 :
  fff = "gpl/m0.001524_Rhox_0.0_outcorr.gpl"
  no_config = 184
  plot_cxorr(fff, no_config, "bo","Charge =0   (nstats = 101)" , 0.1, "RHO_E0") 

if 1 :
  fff = "gpl/m0.001524_Rhox_0.202_outcorr.gpl"
  no_config = 184
  plot_cxorr(fff, no_config, "go","Charge =0.202   (nstats = 101)" , 0.2, "RHO_E") 




##
##  setup plots
##

##plt.title("Pseuodo-scalar meson correlator")
plt.title("Vector meson correlator ")
##plt.legend()



plt.yscale('log')
plt.xlabel('t/a')
plt.ylabel('corr(t)')
plt.xlim(5,15)
##plt.xlim(0,0.25)



##plt.ylim(0,8)
plt.legend()

##plt.savefig("figures/rho_corr.png")
plt.show()

