"""
   Routines to help reformat the correlators

"""

import re




#
#  extract the sweep number from the file name
#  sloppy inversions

def get_config_sloppy(fname,file_tag) :
    cfg = fname.replace(file_tag, "")
    cfg = re.sub("_new_hightol" , "" , cfg)
    cfg = re.sub("_[0-9]+\." , "" , cfg)
    return cfg

#
#  extract the sweep number from the file name
#  precise inversions

def get_config_precise(fname,file_tag) :
    cfg = fname.replace(file_tag, "")
    cfg = re.sub("_newnew_hightol" , "" , cfg)
    cfg = re.sub("_[0-9]+\." , "" , cfg)
    return cfg




#
#  Create a storage area for parametrs
#
#

import sys

class param_choice:
  def __init__(self):
      self.corr_l = dict()
      self.sloppy_file_tag = dict()
      self.prec_file_tag = dict()
      self.o_tag = dict()
      self.dir  = dict()

  def add_option(self,tag,corr_l_, sloppy_file_tag_, prec_file_tag_, o_tag_, dir_) :

      self.corr_l[tag] = corr_l_
      self.sloppy_file_tag[tag] = sloppy_file_tag_
      self.prec_file_tag[tag] = prec_file_tag_
      self.o_tag[tag] =  o_tag_
      self.dir[tag] = dir_


  def loop_option(self) :

     for kk in self.corr_l :
        print ("TAG = ", kk)
        print (self.corr_l[kk])
        print (self.sloppy_file_tag[kk])
        print (self.prec_file_tag[kk])
        print (self.o_tag[kk])
        print (self.dir[kk])


  def get_param(self, tag) :
      
      if not ( tag in self.corr_l   )  :
         print("ERROR tag = " , tag, "not found")
         sys.exit(0)

      return self.corr_l[tag],self.sloppy_file_tag[tag], self.prec_file_tag[tag],  \
             self.o_tag[tag],self.dir[tag]


def setup_param(wdir_, mass_) :


   pc = param_choice()

############################################################

#   pseudo-scalar -- needs to be reformatted

   if 0 :
       corr_l = [  "cc_p_p_d_d_m0.001524_q0.0_m0.001524_q0.0_p000" ]
       sloppy_file_tag = "corrfile.G5G5_charge0.0_t0"
       prec_file_tag = "corrfile.G5G5_charge0.0_fine_t0"
       o_tag = "pseudoscalar_"
   elif 0 :
       corr_l = [  "cc_p_p_d_d_m0.001524_q0.202_m0.001524_q0.202_p000" ]
       sloppy_file_tag = "corrfile.G5G5_charge0.202_t0"
       prec_file_tag = "corrfile.G5G5_charge0.202_fine_t0"
       o_tag = "pseudoscalar_"
   elif 0 :
       corr_l = [  "cc_p_p_d_d_m0.001524_q-0.202_m0.001524_q-0.202_p000" ]
       sloppy_file_tag = "corrfile.G5G5_charge-0.202_t0"
       prec_file_tag = "corrfile.G5G5_charge-0.202_fine_t0"
       o_tag = "pseudoscalar_"


#else:
#  print("Error selection")
#  sys.exit(0)


############################################################

  
##  pseudo-scalar

   ptag = mass_ + "_pion_"

   pc.add_option("UP_G5G5_P","cc_p_p_d_d_m0.001524_q0.202_m0.001524_q0.202_p000",
              "corrfile.G5G5_charge0.202_t0",
              "corrfile.G5G5_charge0.202_fine_t0", ptag, wdir_)

   pc.add_option("UP_G5G5_M","cc_p_p_d_d_m0.001524_q-0.202_m0.001524_q-0.202_p000",
              "corrfile.G5G5_charge-0.202_t0",
              "corrfile.G5G5_charge-0.202_fine_t0", ptag, wdir_)

   pc.add_option("UP_G5G5_0","cc_p_p_d_d_m0.001524_q0.0_m0.001524_q0.0_p000",
              "corrfile.G5G5_charge0.0_t0",
              "corrfile.G5G5_charge0.0_fine_t0", ptag, wdir_)



   pc.add_option("DOWN_G5G5_P","cc_p_p_d_d_m0.003328_q0.101_m0.003328_q0.101_p000",
              "corrfile.G5G5_charge0.101_t0",
              "corrfile.G5G5_charge0.101_fine_t0", ptag, wdir_)

   pc.add_option("DOWN_G5G5_M","cc_p_p_d_d_m0.003328_q-0.101_m0.003328_q-0.101_p000",
              "corrfile.G5G5_charge-0.101_t0",
              "corrfile.G5G5_charge-0.101_fine_t0", ptag, wdir_)

   pc.add_option("DOWN_G5G5_0","cc_p_p_d_d_m0.003328_q0.0_m0.003328_q0.0_p000",
              "corrfile.G5G5_charge0.0_t0",
              "corrfile.G5G5_charge0.0_fine_t0", ptag, wdir_)



   pc.add_option("STRANGE_G5G5_P","cc_p_p_d_d_m0.0677_q0.202_m0.0677_q0.202_p000",
              "corrfile.G5G5_charge0.202_t0",
              "corrfile.G5G5_charge0.202_fine_t0", ptag, wdir_)

   pc.add_option("STRANGE_G5G5_M","cc_p_p_d_d_m0.0677_q-0.202_m0.0677_q-0.202_p000",
              "corrfile.G5G5_charge-0.202_t0",
              "corrfile.G5G5_charge-0.202_fine_t0", ptag, wdir_)

   pc.add_option("STRANGE_G5G5_0","cc_p_p_d_d_m0.0677_q0.0_m0.0677_q0.0_p000",
              "corrfile.G5G5_charge0.0_t0",
              "corrfile.G5G5_charge0.0_fine_t0", ptag, wdir_)



##  vector


#
#     up quark
#
   vtag = mass_ + "_rhox_"

   pc.add_option("UP_GX_P","cc_p_V2_d_d_m0.001524_q0.202_m0.001524_q0.202_p000", 
              "corrfile.GXGX_charge0.202_t0",
              "corrfile.GXGX_charge0.202_fine_t0", vtag, wdir_)

   pc.add_option("UP_GX_M","cc_p_V2_d_d_m0.001524_q-0.202_m0.001524_q-0.202_p000",
              "corrfile.GXGX_charge-0.202_t0", 
              "corrfile.GXGX_charge-0.202_fine_t0",vtag, wdir_)


   pc.add_option("UP_GX_0","cc_p_V2_d_d_m0.001524_q0.0_m0.001524_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)



#
#     down quark
#
   vtag = mass_ + "_rhox_"

   pc.add_option("DOWN_GX_P","cc_p_V2_d_d_m0.003328_q0.101_m0.003328_q0.101_p000", 
              "corrfile.GXGX_charge0.101_t0",
              "corrfile.GXGX_charge0.101_fine_t0", vtag, wdir_)

   pc.add_option("DOWN_GX_M","cc_p_V2_d_d_m0.003328_q-0.101_m0.003328_q-0.101_p000",
              "corrfile.GXGX_charge-0.101_t0", 
              "corrfile.GXGX_charge-0.101_fine_t0",vtag, wdir_)


   pc.add_option("DOWN_GX_0","cc_p_V2_d_d_m0.003328_q0.0_m0.003328_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)



#
#     strange quark
#
   vtag = mass_ + "_rhox_"

   pc.add_option("STRANGE_GX_P","cc_p_V2_d_d_m0.0677_q0.202_m0.0677_q0.202_p000", 
              "corrfile.GXGX_charge0.202_t0",
              "corrfile.GXGX_charge0.202_fine_t0", vtag, wdir_)

   pc.add_option("STRANGE_GX_M","cc_p_V2_d_d_m0.0677_q-0.202_m0.0677_q-0.202_p000",
              "corrfile.GXGX_charge-0.202_t0", 
              "corrfile.GXGX_charge-0.202_fine_t0",vtag, wdir_)


   pc.add_option("STRANGE_GX_0","cc_p_V2_d_d_m0.0677_q0.0_m0.0677_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)

#
#  ml*7  Q= 0.101
#

##   vtag = mass_ + "_rhox_Q0.101"
   vtag = mass_ + "_rhox_"

   pc.add_option("7ML_GX_0.101P","cc_p_V2_d_d_m0.016982_q0.101_m0.016982_q0.101_p000", 
              "corrfile.GXGX_charge0.101_t0",
              "corrfile.GXGX_charge0.101_fine_t0", vtag, wdir_)

   pc.add_option("7ML_GX_0.101M","cc_p_V2_d_d_m0.016982_q-0.101_m0.016982_q-0.101_p000",
              "corrfile.GXGX_charge-0.101_t0", 
              "corrfile.GXGX_charge-0.101_fine_t0",vtag, wdir_)


   pc.add_option("7ML_GX_0","cc_p_V2_d_d_m0.016982_q0.0_m0.016982_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)



#
#  ml*7  Q= 0.202
#

##   vtag = mass_ + "_rhox_Q0.202"
   vtag = mass_ + "_rhox_"

   pc.add_option("7ML_GX_0.202P","cc_p_V2_d_d_m0.016982_q0.202_m0.016982_q0.202_p000", 
              "corrfile.GXGX_charge0.202_t0",
              "corrfile.GXGX_charge0.202_fine_t0", vtag, wdir_)

   pc.add_option("7ML_GX_0.202M","cc_p_V2_d_d_m0.016982_q-0.202_m0.016982_q-0.202_p000",
              "corrfile.GXGX_charge-0.202_t0", 
              "corrfile.GXGX_charge-0.202_fine_t0",vtag, wdir_)



#
#  5*ml  Q=0.101
#

##   vtag = mass_ + "_rhox_Q0.101"
   vtag = mass_ + "_rhox_"

   pc.add_option("5ML_GX_0.101P","cc_p_V2_d_d_m0.01213_q0.101_m0.01213_q0.101_p000", 
              "corrfile.GXGX_charge0.101_t0",
              "corrfile.GXGX_charge0.101_fine_t0", vtag, wdir_)

   pc.add_option("5ML_GX_0.101M","cc_p_V2_d_d_m0.01213_q-0.101_m0.01213_q-0.101_p000",
              "corrfile.GXGX_charge-0.101_t0", 
              "corrfile.GXGX_charge-0.101_fine_t0",vtag, wdir_)


   pc.add_option("5ML_GX_0","cc_p_V2_d_d_m0.01213_q0.0_m0.01213_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)

#
#  5*ml  Q=0.202
#

##   vtag = mass_ + "_rhox_Q0.202"
   vtag = mass_ + "_rhox_"

   pc.add_option("5ML_GX_0.202P","cc_p_V2_d_d_m0.01213_q0.202_m0.01213_q0.202_p000", 
              "corrfile.GXGX_charge0.202_t0",
              "corrfile.GXGX_charge0.202_fine_t0", vtag, wdir_)

   pc.add_option("5ML_GX_0.202M","cc_p_V2_d_d_m0.01213_q-0.202_m0.01213_q-0.202_p000",
              "corrfile.GXGX_charge-0.202_t0", 
              "corrfile.GXGX_charge-0.202_fine_t0",vtag, wdir_)




#
#  3*ml  Q= 0.101
#

##   vtag = mass_ + "_rhox_Q0.101"
   vtag = mass_ + "_rhox_"

   pc.add_option("3ML_GX_0.101P","cc_p_V2_d_d_m0.007278_q0.101_m0.007278_q0.101_p000", 
              "corrfile.GXGX_charge0.101_t0",
              "corrfile.GXGX_charge0.101_fine_t0", vtag, wdir_)

   pc.add_option("3ML_GX_0.101M","cc_p_V2_d_d_m0.007278_q-0.101_m0.007278_q-0.101_p000",
              "corrfile.GXGX_charge-0.101_t0", 
              "corrfile.GXGX_charge-0.101_fine_t0",vtag, wdir_)


   pc.add_option("3ML_GX_0","cc_p_V2_d_d_m0.007278_q0.0_m0.007278_q0.0_p000",
               "corrfile.GXGX_charge0.0_t0",
                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)



#
#  3*ml  Q= 0.202
#

##   vtag = mass_ + "_rhox_Q0.202"
   vtag = mass_ + "_rhox_"

   pc.add_option("3ML_GX_0.202P","cc_p_V2_d_d_m0.007278_q0.202_m0.007278_q0.202_p000", 
              "corrfile.GXGX_charge0.202_t0",
              "corrfile.GXGX_charge0.202_fine_t0", vtag, wdir_)

   pc.add_option("3ML_GX_0.202M","cc_p_V2_d_d_m0.007278_q-0.202_m0.007278_q-0.202_p000",
              "corrfile.GXGX_charge-0.202_t0", 
              "corrfile.GXGX_charge-0.202_fine_t0",vtag, wdir_)


#   pc.add_option("3ML_GX_0","cc_p_V2_d_d_m0.007278_q0.0_m0.007278_q0.0_p000",
#               "corrfile.GXGX_charge0.0_t0",
#                "corrfile.GXGX_charge0.0_fine_t0", vtag , wdir_)





   pc.loop_option()

   return pc
