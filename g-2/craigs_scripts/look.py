import os
import shutil

def get_config_list_tag(gdir, tag) :

    clist = [] 
#    print ('Searching for configurations in ' , gdir)
    for filename in os.listdir(gdir):
        if tag in filename :
           clist.append(filename)

    slist = sorted(clist)
    return slist

#
#
#

def move_bad(gdir, gdir_bad, cfg_list) :

   
    for cfg_ in cfg_list :
       for filename in os.listdir(gdir):
          if str(cfg_) in filename :
            ff_ = gdir + "/" + filename
            ff_out = gdir_bad + "/" + filename
            print("Move " , ff_)
            shutil.move(ff_, ff_out)

##
##
##


def get_cfg_list(dir_, tag_) :

  dd = get_config_list_tag(dir_ , tag_ ) 

  cfg_ = [] 
  for d_ in dd :
     tmp = d_.split(".")
     tmp_ = tmp[-1].replace("_new_hightol", "")
     tmp_ = tmp_.replace("_newnew_hightol", "")
     cfg_.append(int(tmp_))

  cfg_list = sorted(list(set(cfg_)))

  print("Tag " , tag_, "dir" , dir_)
  print("Files found " , len(dd))
  print("Number of configurations = " , len(cfg_list))
  return cfg_list

##mass = "m0.001524"
##mass = "m0.0677"
##charge="0.202"

##charge="0.101"
##mass = "m0.003328"

##  -----  3*ml
charge="0.101"
mass = "m0.007278"

##charge="0.202"
##mass = "m0.007278"


##  -----  5*ml
charge="0.101"
mass = "m0.01213"

##charge="0.202"
##mass = "m0.01213"

##  -----  7*ml
#charge="0.101"
#mass = "m0.01698"

##charge="0.202"
##mass = "m0.01698"





dir = "gdata/" + mass
tag = "corrfile.GXGX_charge-" + charge +   "_t0"

dir_bad = "gdata/" + mass + "_bad"


cfg_listA = get_cfg_list(dir, "corrfile.GXGX_charge-"+ charge+"_t0")
cfg_listB = get_cfg_list(dir, "corrfile.GXGX_charge"+ charge+"_t0")
cfg_listC = get_cfg_list(dir, "corrfile.GXGX_charge0.0_t0")

cfg_listA_f = get_cfg_list(dir, "corrfile.GXGX_charge-"+ charge+"_fine_t0")
cfg_listB_f = get_cfg_list(dir, "corrfile.GXGX_charge"+ charge+"_fine_t0")
cfg_listC_f = get_cfg_list(dir, "corrfile.GXGX_charge0.0_fine_t0")


cfg_list_ = set(cfg_listA).intersection(set(cfg_listB)) 
cfg_list_ = set(cfg_listA).intersection(set(cfg_listC_f)) 
cfg_list_ = set(cfg_listA).intersection(set(cfg_listB_f)) 
cfg_list_ = set(cfg_listA).intersection(set(cfg_listA_f)) 
cfg_list_ = set(cfg_list_).intersection(set(cfg_listC)) 

cfg_list = list(cfg_list_) 


cfg_list_t = cfg_listA + cfg_listB + cfg_listC + cfg_listA_f  +  cfg_listB_f  + cfg_listC_f 
cfg_list_t = list(set(cfg_list_t))

print("Total number of configurations " , len(cfg_list_t))
print("Number of common configurations " , len(cfg_list))

cfg_missing = list((set(cfg_list_t).difference(cfg_list)) )

print("Missing values in first list:", cfg_missing)

for cfg in cfg_missing :
   print(cfg)


move_bad(dir, dir_bad, cfg_missing) 


