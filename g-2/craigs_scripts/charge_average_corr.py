#!/Users/cmcneile/anaconda3/bin/python3.6

import pandas  as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("/Users/cmcneile/projects/qed/pythonlib")
import store_results
import pcorr

from param import *
##from analysis import *

##
##  load the two files
##

def charge_average(name_av, corr_tag_P, name_P, corr_tag_M, name_M) :
  
  try:
    df_sloppy_P, qmass = pcorr.load_corr(name_P,corr_tag_P)
  except:
    print("Error loading " , name_P)
    sys.exit(0)

  print(df_sloppy_P.head())

  print("--------------------------------------------------")
  try:
    df_sloppy_M, qmass = pcorr.load_corr(name_M,corr_tag_M)
  except:
    print("Error loading " , name_M)
    sys.exit(0)

  print(df_sloppy_M.head())

  ##
  print("Averaging the dataframes")

  df_sloppy_M = df_sloppy_M.add(df_sloppy_P)
  df_sloppy_M /= 2
  df_sloppy_M['sweep'] = df_sloppy_M['sweep'].astype(int)

  print(df_sloppy_M.head())

  print("Writing averaged sloppy correlators to " , name_av )
  df_sloppy_M.to_csv(name_av)



##################################################
Q = "0.202"
##Q = "0.101"

if 1 : 
   corr_tag_P = "cc_p_p_d_d_" + mass + "_q" +Q+  "_" + mass + "_q" + Q + "_p000"
   corr_tag_M = "cc_p_p_d_d_" + mass + "_q-" + Q + "_" + mass +  "_q-" + Q + "_p000"

# ----  sloppy ----
   name_av = "data/" + mass + "_pion_corrfile.G5G5_charge_AV_" + Q + "_t0.csv"
   name_P = "data/" + mass + "_pion_corrfile.G5G5_charge" + Q + "_t0.csv"
   name_M = "data/" + mass + "_pion_corrfile.G5G5_charge-" + Q + "_t0.csv"

   name_av_prec = "data/" + mass + "_pion_corrfile.G5G5_charge_AV_" + Q + "_fine_t0.csv"
   name_P_prec = "data/" + mass + "_pion_corrfile.G5G5_charge" + Q + "_fine_t0.csv"
   name_M_prec = "data/" + mass + "_pion_corrfile.G5G5_charge-" + Q + "_fine_t0.csv"

elif 0 : 
   corr_tag_P = "cc_p_V2_d_d_" + mass + "_q" + Q + "_" + mass + "_q" + Q + "_p000"
   corr_tag_M = "cc_p_V2_d_d_" + mass + "_q-" + Q + "_" + mass + "_q-" + Q + "_p000"

##   mass = "m0.001524"

# ----  sloppy ----
   name_av = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_" + Q + "_t0.csv"
   name_P  = "data/" + mass + "_rhox_corrfile.GXGX_charge" + Q + "_t0.csv"
   name_M  = "data/" + mass + "_rhox_corrfile.GXGX_charge-" + Q + "_t0.csv"

   name_av_prec = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_" + Q + "_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge" + Q + "_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge-" + Q + "_fine_t0.csv"



elif 0: 
#      ----  up ----

##   mass = "m0.007278"

   corr_tag_P = "cc_p_V2_d_d_" + mass + "_q0.101_" + mass + "_q0.101_p000"
   corr_tag_M = "cc_p_V2_d_d_" + mass + "_q-0.101_" + mass + "_q-0.101_p000"


# ----  sloppy ----
   name_av = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_0.101_t0.csv"
   name_P  = "data/" + mass + "_rhox_corrfile.GXGX_charge0.101_t0.csv"
   name_M  = "data/" + mass + "_rhox_corrfile.GXGX_charge-0.101_t0.csv"

   name_av_prec = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_0.101_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge0.101_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge-0.101_fine_t0.csv"



elif 0 : 
#      ----  3*ml ----

#   Q="0.101"
   Q="0.202"
   mass = "m0.007278"
   vtag = "rhox_Q" + Q

   corr_tag_P = "cc_p_V2_d_d_m0.007278_q"  + Q +  "_m0.007278_q" + Q +  "_p000"
   corr_tag_M = "cc_p_V2_d_d_m0.007278_q-" + Q + "_m0.007278_q-" + Q +  "_p000"


# ----  sloppy ----
   name_av = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" + Q +  "_t0.csv"
   name_P  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q +  "_t0.csv"
   name_M  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_t0.csv"

   name_av_prec = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" +Q+  "_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q +  "_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_fine_t0.csv"


elif 0 : 
#      ----  5*ml ----

#   Q="0.101"
   Q="0.202"
   mass = "m0.01213"
   vtag = "rhox_Q" + Q

   corr_tag_P = "cc_p_V2_d_d_m0.01213_q" + Q + "_m0.01213_q" + Q + "_p000"
   corr_tag_M = "cc_p_V2_d_d_m0.01213_q-" + Q + "_m0.01213_q-" + Q + "_p000"


# ----  sloppy ----
   name_av = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" + Q + "_t0.csv"
   name_P  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q + "_t0.csv"
   name_M  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_t0.csv"

   name_av_prec = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" + Q + "_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q + "_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_fine_t0.csv"



elif 0 : 
#      ----  7*ml ----

##   Q="0.101"
   Q="0.202"
   mass = "m0.01698"
##   vtag = "rhox_Q" + Q
   vtag = "rhox_"

   corr_tag_P = "cc_p_V2_d_d_m0.016982_q" + Q + "_m0.016982_q" + Q + "_p000"
   corr_tag_M = "cc_p_V2_d_d_m0.016982_q-" + Q + "_m0.016982_q-" + Q + "_p000"


# ----  sloppy ----
   name_av = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" + Q + "_t0.csv"
   name_P  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q + "_t0.csv"
   name_M  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_t0.csv"

   name_av_prec = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge_AV_" + Q + "_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge" + Q + "_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_" +  vtag + "corrfile.GXGX_charge-" + Q + "_fine_t0.csv"



elif 0 : 
#      ----  strange ----

   mass = "m0.0677"

   corr_tag_P = "cc_p_V2_d_d_" + mass + "_q0.202_" + mass + "_q0.202_p000"
   corr_tag_M = "cc_p_V2_d_d_" + mass + "_q-0.202_" + mass + "_q-0.202_p000"



# ----  sloppy ----
   name_av = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_0.202_t0.csv"
   name_P  = "data/" + mass + "_rhox_corrfile.GXGX_charge0.202_t0.csv"
   name_M  = "data/" + mass + "_rhox_corrfile.GXGX_charge-0.202_t0.csv"

   name_av_prec = "data/" + mass + "_rhox_corrfile.GXGX_charge_AV_0.202_fine_t0.csv"
   name_P_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge0.202_fine_t0.csv"
   name_M_prec  = "data/" + mass + "_rhox_corrfile.GXGX_charge-0.202_fine_t0.csv"






print("--------------------------------------------------")
print("Average +/- correlators")
print("--------------------------------------------------")

charge_average(name_av, corr_tag_P, name_P, corr_tag_M, name_M) 
charge_average(name_av_prec, corr_tag_P, name_P_prec, corr_tag_M, name_M_prec) 

sys.exit(0)
