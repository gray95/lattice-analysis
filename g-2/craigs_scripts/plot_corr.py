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
from analysis import *

##################################################

try:
  df_sloppy, qmass = pcorr.load_corr(name,corr_tag)
except:
  print("Error loading " , name)
  sys.exit(0)

df_sloppy_store = df_sloppy.copy()

#
#  precise time sources
#

try:
  df_p, qmass_p = pcorr.load_corr(name_p,corr_tag)
except:
  print("Error loading "  , name_p)
  sys.exit(0)

#
#  check for missing sweeps
#


XX_p = df_p['sweep']
XX = df_sloppy['sweep']

print("Missing values in second list:", (set(XX.values).difference(XX_p.values))) 
print("Missing values in first list:", (set(XX_p.values).difference(XX.values))) 



##################################################
#  Extract the t_src at precise

tsrc_prec = df_p[[ 'tsrc' , 'sweep' ] ]
print(tsrc_prec)
tsrc_prec.set_index('sweep', inplace=True)
ggg = tsrc_prec.to_dict()
print(ggg)
t_prec = ggg['tsrc']
##print(t_prec)


##################################################

# create sloppy at same time source as precise
#
df_sloppy_at_prec = df_sloppy.copy()

drop_store_notprec = [] 
drop_store_at_prec = [] 

for index, row in df_sloppy_at_prec.iterrows():
#    print(row['sweep'], row['tsrc'])
#    print(type(row['sweep']), type(row['tsrc']))

    sweep_ = int(row['sweep'])
    tsrc_  = int(row['tsrc'])
    tsrc_pp_ = t_prec[sweep_]

    
    if tsrc_ == tsrc_pp_  :
       print("match" , tsrc_ , "index = " , index)
       drop_store_at_prec.append(index)
    else:
       drop_store_notprec.append(index)
       print("NOT MATCH" , tsrc_ , "index = " , index)

##sys.exit(0)

#
#  sloppy at tsrc of precise 
#
df_sloppy_at_prec.drop(df_sloppy_at_prec.index[drop_store_notprec] , inplace=True )
df_sloppy_at_prec.reset_index(drop=True,  inplace=True)
df_sloppy_at_prec.drop(["tsrc"], inplace=True, axis=1)
df_sloppy_at_prec.set_index('sweep' , inplace=True )

##df_sloppy_at_prec.drop(["Unnamed: 0"], inplace=True, axis=1)

# 
# remove tsrc at the precise srcf for sloppy
#

df_sloppy.drop(df_sloppy.index[drop_store_at_prec] , inplace=True )
df_sloppy.reset_index(drop=True,  inplace=True)
df_sloppy.drop(["tsrc"], inplace=True, axis=1)

df_sloppy.set_index('sweep')
# average over the sources
df_sloppy_av  = df_sloppy.groupby("sweep").mean()
##df_sloppy_av.drop(["Unnamed: 0"], inplace=True, axis=1)  
df_sloppy_av.sort_index(inplace=True, ascending=False)
##df_sloppy_av.reset_index(drop=True,  inplace=True)

##
##  process the precise data frame
##

df_p.drop(["tsrc"], inplace=True, axis=1)
df_p.reset_index(drop=True,  inplace=True)
df_p.set_index('sweep', inplace=True)
##df_p.drop(["Unnamed: 0"], inplace=True, axis=1)  



##################################################
print("df_sloppy_av (no precise)")
print(df_sloppy_av.head())

print("df_sloppy_at_prec")
print(df_sloppy_at_prec.head())

print("df_p")
print(df_p.head())



############################################################
#   craete the corrected 
############################################################

df_corrected = df_p.copy()
df_corrected += df_sloppy_av  -  df_sloppy_at_prec
##df_corrected = df_sloppy_av.copy()  ## debug
##df_corrected = df_sloppy_at_prec.copy()  ## debug

##df_corrected = df_sloppy_av.copy()

##sys.exit(0)

#  debug debug debug debug 
##df_corrected = df_sloppy_av.copy()
#df_corrected = df_sloppy_store.copy()
#df_corrected = df_corrected.groupby("sweep").mean()
#  debug debug debug debug 

#
print("df_corrected")
print(df_corrected.head() )

##sys.exit(0)

##print("Precise time sources = " , df_p["tsrc"].unique())

##df_p.drop(["tsrc"], inplace=True, axis=1)
#df_p.drop(["corr_tag"], inplace=True, axis=1)
#df_p.drop(["qmass"], inplace=True, axis=1)
##df_p.set_index('sweep' , inplace=True )

if gpl_tag :

   nrm = 3.0*32*32*32
#   nrm = 1
   df_corrected /= nrm
   pcorr.dump_corr_lepage_format(df_corrected,nt, "gpl/" + mass + "_" + ptitle + "_outcorr.gpl", gpl_tag) 
#   df_corrected.to_csv(mass + "_" + ptitle + "_corrected.csv")

   #  write in csv format
   #   df_p /= nrm
   # print(df_p.head())
   # pcorr.dump_corr_lepage_format(df_p,nt, "gpl/prec_" + mass + "_" + ptitle + "_outcorr.gpl", gpl_tag) 

##sys.exit(0)



##  --------------------------------------------------

ttt       = dict()
corr      = dict()
corr_err  = dict()

ttt["t0"],corr["t0"],corr_err["t0"] = pcorr.get_mean_corr(df_sloppy_at_prec,nt,0.05)
ttt["av"],corr["av"],corr_err["av"] = pcorr.get_mean_corr(df_sloppy_av,nt)

ttt["t_prec"],corr["t_prec"],corr_err["t_prec"] = pcorr.get_mean_corr(df_p,nt,-0.05)

ttt["t_corrected"],corr["t_corrected"],corr_err["t_corrected"] = pcorr.get_mean_corr(df_corrected,nt,-0.1)


#
#
#
noconfig = len(XX_p.values)
print("Qmass = " , qmass)
print("Number of configs = " , noconfig)


store_results.to_pickle("pickle/" + ptitle + "_corr.pickle", ttt,corr, corr_err, qmass, noconfig)



##
##  plot
##

do_corr_plot = lambda tag__, mk, lab : plt.errorbar(ttt[tag__],corr[tag__],corr_err[tag__], fmt=mk , label=lab)


def do_ratio(tag__, mk, lab) :
  rat = np.array(corr_err[tag__]) /  np.array( corr[tag__])
  plt.plot(ttt[tag__],rat, mk , label=lab)



##do_plot = do_corr_plot 
do_plot = do_ratio


do_plot("t0", "ro", "(sloppy) single tsrc") 
do_plot("t_prec", "bo", "(precise) single tsrc") 
do_plot("av", "go", "(sloppy) averaged tsrc") 
do_plot("t_corrected", "yo", "corrected") 

##plt.title(ptitle + " correlator at quark mass = " + str(qmass) )
plt.title(ptitle + " error/mean correlator " )
plt.legend()

plt.xlabel('t/a')
plt.ylabel(r'$\sigma(t)$/corr(t)')
##plt.ylabel(r'corr(t)')
##plt.xlim(2.8,3.8)
plt.xlim(1.8,5.8)

plt.yscale('log')
##plt.ylim(100,1400)
##plt.ylim(0,500)


plt.savefig("figures/" + ptitle + "_corr.png")

plt.show()
