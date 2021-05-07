#!/usr/local/bin/python3


import math



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


def plot_cxorr(fff_, no_config_, fmt_, lab_, sh) :
  corr_tag = "RHO_E_AV"
  corr = util.load_data(fff, nt,no_config_, corr_tag)

  print("Computing jackknife correlators")
  tt, corr_mean, corr_err = util.calc_corr(corr, nt,no_config_,sh ) 
  plt.errorbar(tt, corr_mean , corr_err  ,  fmt= fmt_ , label=lab_ )


#
#  
#

if 1 :
  fff = "gpl/m0.001524_Rhox_0.202_chargeAV_outcorr.gpl"
  no_config = 101
  plot_cxorr(fff, no_config, "ro","up  (101 lattices)" , 0.0) 


if 1 :
  fff = "gpl/m0.0677_Rhox_0.202_chargeAV_outcorr.gpl"
  no_config = 101
  plot_cxorr(fff, no_config, "go","strange  (101 lattices)" , 0.1) 

fff = "gpl/m0.01698_Rhox_0.101_chargeAV_outcorr.gpl"
no_config = 121
plot_cxorr(fff, no_config, "bo","7*ml  (121 lattices)" , 0.2) 



##
##  setup plots
##

##plt.title("Pseuodo-scalar meson correlator")
plt.title("Neutral vector meson correlator (Q=0.202)")
##plt.legend()



plt.yscale('log')
plt.xlabel('t/a')
plt.ylabel('corr(t)')
plt.xlim(0,20)
##plt.xlim(0,0.25)



##plt.ylim(0,8)
plt.legend()

plt.savefig("figures/ml_rho_corr.png")
plt.show()

