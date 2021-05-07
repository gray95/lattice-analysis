"""
   Analysis choices
"""

import sys

gpl_tag = ""

# perhaps pull from param.h
mass = "m0.001524"
##mass = "m0.003328"
##mass = "m0.0677"

##mass = "m0.01698"
##mass = "m0.01213"
##mass = "m0.007278"

Q = "0.202"
##Q = "0.101"

if 0 :
  name = "data/" +  mass + "_pion_corrfile.G5G5_charge0.0_t0.csv"
  name_p = "data/" +  mass + "_pion_corrfile.G5G5_charge0.0_fine_t0.csv"
  corr_tag = "cc_p_p_d_d_" + mass + "_q0.0_" + mass + "_q0.0_p000"
  ptitle = "Pseuodoscalar_0.0"
  gpl_tag = "PSEUDO_E0"
elif 0 :
  name = "data/" +  mass + "_pion_corrfile.G5G5_charge" + Q + "_t0.csv"
  name_p = "data/" +  mass + "_pion_corrfile.G5G5_charge" + Q + "_fine_t0.csv"
  corr_tag = "cc_p_p_d_d_" + mass + "_q" + Q + "_" + mass + "_q" + Q + "_p000"
  ptitle = "Pseuodoscalar_" + Q
  gpl_tag = "PSEUDO_EP"
elif 0:
  name = "data/" +  mass + "_pion_corrfile.G5G5_charge-0.202_t0.csv"
  name_p = "data/" +  mass + "_pion_corrfile.G5G5_charge-0.202_fine_t0.csv"
  corr_tag = "cc_p_p_d_d_m0.001524_q-0.202_m0.001524_q-0.202_p000"
  ptitle = "Pseuodoscalar_-0.202"
  gpl_tag = "PSEUDO_EM"
elif 0 :
  name = "data/" +  mass + "_pion_corrfile.G5G5_charge_AV_" + Q + "_t0.csv"
  name_p = "data/" +  mass + "_pion_corrfile.G5G5_charge_AV_" + Q + "_fine_t0.csv"
  corr_tag = ""
  ptitle = "Pseuodoscalar_" + Q  + "_chargeAV"
  gpl_tag = "PSEUDO_E_AV"
#
#    Vector
#
elif 1 :
  name = "data/" +  mass + "_rhox_corrfile.GXGX_charge_AV_" + Q + "_t0.csv"
  name_p = "data/" +  mass + "_rhox_corrfile.GXGX_charge_AV_" + Q + "_fine_t0.csv"
  corr_tag = ""
  ptitle = "Rhox_" + Q + "_chargeAV"
  gpl_tag = "RHO_E_AV"

elif 0 :
  name = "data/" +  mass + "_rhox_corrfile.GXGX_charge" + Q + "_t0.csv"
  name_p = "data/" +  mass + "_rhox_corrfile.GXGX_charge" + Q + "_fine_t0.csv"
  corr_tag = "cc_p_V2_d_d_" + mass + "_q" + Q + "_" + mass + "_q" + Q + "_p000"
  ptitle = "Rhox_" + Q
  gpl_tag = "RHO_EP"
elif 0 :
  name = "data/" +mass + "_rhox_corrfile.GXGX_charge0.0_t0.csv"
  name_p = "data/" + mass + "_rhox_corrfile.GXGX_charge0.0_fine_t0.csv"
  corr_tag = "cc_p_V2_d_d_" + mass + "_q0.0_" + mass + "_q0.0_p000"
  ptitle = "Rhox_0.0"
  gpl_tag = "RHO_E0"


elif 0 :
##  7Ml 
##   additional digit in tag

  name = "data/" +mass + "_rhox_corrfile.GXGX_charge0.0_t0.csv"
  name_p = "data/" + mass + "_rhox_corrfile.GXGX_charge0.0_fine_t0.csv"
  corr_tag = "cc_p_V2_d_d_m0.016982_q0.0_m0.016982_q0.0_p000"
  ptitle = "Rhox_0.0"
  gpl_tag = "RHO_E0"


#
#        down
#
elif 0 :
  name = "data/" +  mass + "_rhox_corrfile.GXGX_charge_AV_0.101_t0.csv"
  name_p = "data/" +  mass + "_rhox_corrfile.GXGX_charge_AV_0.101_fine_t0.csv"
  corr_tag = ""
  ptitle = "Rhox_0.101_chargeAV"
  gpl_tag = "RHO_E_AV"
##
elif 0 :
  name = "data/" +mass + "_rhox_Q0.101corrfile.GXGX_charge0.0_t0.csv"
  name_p = "data/" + mass + "_rhox_Q0.101corrfile.GXGX_charge0.0_fine_t0.csv"
  corr_tag = ""
  ptitle = "Rhox_0.0"
  gpl_tag = "RHO_E0"

elif 0 :
  name = "data/" +mass + "_rhox_Q0.101corrfile.GXGX_charge_AV_0.101_t0.csv"
  name_p = "data/" + mass + "_rhox_Q0.101corrfile.GXGX_charge_AV_0.101_fine_t0.csv"
  corr_tag = ""
  ptitle = "Rhox_0.101_chargeAV"
  gpl_tag = "RHO_E_AV"

elif 0 :
  name = "data/" +mass + "_rhox_Q0.202corrfile.GXGX_charge_AV_0.202_t0.csv"
  name_p = "data/" + mass + "_rhox_Q0.202corrfile.GXGX_charge_AV_0.202_fine_t0.csv"
  corr_tag = ""
  ptitle = "Rhox_0.202_chargeAV"
  gpl_tag = "RHO_E_AV"

try: 
  print("Corr tag = " , corr_tag)
except:
  print("Error corr_tag not defined")
  sys.exit(0)
