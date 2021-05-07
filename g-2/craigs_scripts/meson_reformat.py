#!/Users/cmcneile/anaconda3/bin/python3.6

import os
import re
import sys

sys.path.append("/Users/cmcneile/projects/qed/pythonlib")
import util
from param import *
from setup_reformat import *

############################################################

##
##  Vector
##

pc = setup_param(wdir_fine, mass)


tag_list = dict()

#tag_list["m0.001524"] = [ "UP_G5G5_P" , "UP_G5G5_M" , "UP_G5G5_0"  ]
#tag_list["m0.003328"] = [ "DOWN_G5G5_0" , "DOWN_G5G5_P" , "DOWN_G5G5_M"  ]
#tag_list["m0.0677"] = [ "STRANGE_G5G5_0" , "STRANGE_G5G5_P" , "STRANGE_G5G5_M"  ]

tag_list["m0.001524"] = [ "UP_G5G5_P" , "UP_G5G5_M" , "UP_G5G5_0" , "UP_GX_0" , "UP_GX_P" , "UP_GX_M" ]  
tag_list["m0.003328"] = [ "DOWN_G5G5_0" , "DOWN_G5G5_P" , "DOWN_G5G5_M" , "DOWN_GX_0" , "DOWN_GX_P" , "DOWN_GX_M"   ]

tag_list["m0.0677"] = [ "STRANGE_G5G5_0" , "STRANGE_G5G5_P" , "STRANGE_G5G5_M" , "STRANGE_GX_0" , "STRANGE_GX_P" , "STRANGE_GX_M" ]

tag_list["m0.007278"] = [ "3ML_GX_0.101P" , "3ML_GX_0.101M" , "3ML_GX_0" , "3ML_GX_0.202P" , "3ML_GX_0.202M"  ]
tag_list["m0.01213"] = [ "5ML_GX_0.101P" , "5ML_GX_0.101M" , "5ML_GX_0" , "5ML_GX_0.202P" , "5ML_GX_0.202M" ]
tag_list["m0.01698"] = [ "7ML_GX_0.101P" , "7ML_GX_0.101M" , "7ML_GX_0" , "7ML_GX_0.202P" , "7ML_GX_0.202M" ]   


##  --------------------------------------------------
##tag_list["m0.001524"] = [ "UP_GX_0" , "UP_GX_P" , "UP_GX_M"  ]
##tag_list["m0.0677"] = [ "STRANGE_GX_0" , "STRANGE_GX_P" , "STRANGE_GX_M"  ]
##tag_list["m0.003328"] = [ "DOWN_GX_0" , "DOWN_GX_P" , "DOWN_GX_M"  ]

##tag_list["m0.01698"] = [ "7ML_GX_0.101P" , "7ML_GX_0.101M" , "7ML_GX_0" ]
##tag_list["m0.01698"] = [ "7ML_GX_0.202P" , "7ML_GX_0.202M" ]

#tag_list["m0.01213"] = [ "5ML_GX_0.101P" , "5ML_GX_0.101M" , "5ML_GX_0" ]
##tag_list["m0.01213"] = [ "5ML_GX_0.202P" , "5ML_GX_0.202M"  ]

##tag_list["m0.007278"] = [ "3ML_GX_0.101P" , "3ML_GX_0.101M" , "3ML_GX_0" ]


for tag in  tag_list[mass] :
 
   print("------------------------------------------------------------")
   print("Working on " , tag )
   print("------------------------------------------------------------")

   corr_l,sloppy_file_tag,prec_file_tag,o_tag,wdir = pc.get_param(tag)


############################################################
##  Precise correlators
############################################################

   print("Extracting the precise meson correlators")

   cfglist_precise =  util.get_config_list_tag(wdir, prec_file_tag) 
   print(cfglist_precise)
   print("INFO Number of fine correlators ",  len(cfglist_precise) )

   fout = data_dir + "/" + o_tag + prec_file_tag + ".csv"
   util.write_corr_from_files_csv(cfglist_precise, get_config_precise, wdir, corr_l, prec_file_tag, nt, fout, False ) 

##sys.exit(0)

############################################################
##  Sloppy correlators
############################################################

   print("Extracting the sloppy correlators")

   cfglist_coarse =  util.get_config_list_tag(wdir, sloppy_file_tag) 

   cfg_list = [] 
   for cfg in cfglist_coarse :
     cfg_ =  get_config_sloppy(cfg,sloppy_file_tag)
     cfg_list.append(cfg_)

   tmp = set(cfg_list)
   print("Number of lattices with sloppy measurements  = " , len(tmp))

##sys.exit(0)

   fout = data_dir + "/" + o_tag + sloppy_file_tag + ".csv"
   util.write_corr_from_files_csv(cfglist_coarse,get_config_sloppy, wdir, corr_l,sloppy_file_tag,nt,fout, False ) 


