import os
import re
import sys

#
#  Read hadron data from a MILC data file
#
def extract_data(filename_in,tag):

#  print ("Reading data from  " , filename_in)
  try:
    f = open(filename_in, 'r')
  except:
    print ("Error opening " , filename_in)
    sys.exit(1)

  store = [] 

  ct = 0 
  dump = 0  

  for line in f:
    ct = ct + 1
    
    ll  = line.rstrip('\n')

    if re.search( r'---'  ,   ll ) :
      dump = 0
      ct = 0 

    if dump  and ct > 2 :
      tmp = ll.split()
#      print(tmp[1])
      store.append(float(tmp[1]  ) )

    if re.search( tag  ,   ll ) :
      dump =  1
      ct = 0 

    ct = ct + 1

  f.close



  return store


def get_config_list(gdir) :

    clist = [] 
    print ('Searching for configurations in ' , gdir)
    for filename in os.listdir(gdir):
           cfg = filename
           clist.append(cfg)

    slist = sorted(clist)
    #print(slist)
    return slist

#################################

def timesym_corr(corr) :

   corrAv  = [] 
   tend    = len(corr)
   corrAv  = corr
   tmid = int( tend / 2 )

   for t in range(1 , tmid ) :
     tneg = tend - t 
     cc = ( corr[t] + corr[tneg] ) * 0.5 
#     print (t , tneg , corr[t] , corr[tneg] , " = " , cc )
     corrAv[t]    = cc
     corrAv[tneg] = cc

##   sys.exit(0) 
   return corrAv

#################################
#	PROGRAM STARTS HERE
#################################

a = 'PION_PS_G_H_c_d_m0.440_m0.440_p000'
#b = 'PION_PS_Hy_A_Gonemp_LG_m0.450_m0.450_p000'
#c = 'PION_PS_Ga_Ga_Hy_GL_m0.450_m0.450_p000'
#d = 'PION_PS_Ga_Ga_onemp_GG_m0.450_m0.450_p000'
corr_tag =  [a]
tag = ["onemp.gg"]

# l3296f211b630m0074m037m440-coul-v5\ks_spectrum_vary
# ks_spectrum_vary_onemp_lrho\corr_t0\200\meson_hybrid_corr.txt
# C:\Users\gsinharay\OneDrive - University of Plymouth\milc data 2019\raw_data\l3296f211b630m0074m037m440-coul-v5\onemp_smearfirst\corr
wdir = os.path.join('raw_data', 'l3296f211b630m0074m037m440-coul-v5', 'onemp_smearfirst', 'corr') 
cfglist =  get_config_list(wdir) 

target = a[:-12] + '.txt' # output file

l = open(target, "w+")
for i in range(len(tag)):   
    for cfg in cfglist:
        fff = os.path.join(wdir, cfg, 'meson_hybrid_corr.txt')
        corr = extract_data(fff,corr_tag[i])
        l.write(tag[i])
        l.write(" ")
        for cc in corr:      
            l.write(("%s" % cc))
            l.write(" ")
        l.write("\n")
  


  
