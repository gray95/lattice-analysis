import numpy as np
import gvar as gv
import sys
import os
sys.path.append('..')
from commonweal import ainv, hbarc

e_str = ['vc', 'c', 'c-10', 'f', 'f-5']

resdir = '/home/gray/Desktop/lattice-analysis/code/out'

res = { x:gv.load(os.path.join(resdir,x,'onemp.p')) for x in e_str } 
onemp_am = { y:res[y]['fit.p']['1mp.dE'][0] for y in res }
onemp_m =  { y:ainv[y]*onemp_am[y] for y in onemp_am }
pp_am = { y:res[y]['fit.p']['pp.dE'][0] for y in res }
pp_m =  { y:ainv[y]*pp_am[y] for y in pp_am }

ratio = { e:onemp_am[e]/pp_am[e] for e in e_str }

E = 'f-5'
print(E)
print("onemp lat",onemp_am[E])
print("onemp GeV",onemp_m[E])
print("PP lat", pp_am[E])
print("PP GeV",pp_m[E])
print("ratio %s"%(onemp_am[E]/pp_am[E]))

#print(gv.corr(onemp_m['f-5'],pp_m['f-5']))


## other groups + old + 1-- paraphernalia
onemp_a = {"hadspec": 0.12, "ma": 0.18, "regensburg": 0.1145}

onemp_mass = {"hadspec": gv.gvar('4.310(29)'),
                                                "ma": gv.gvar('4.309(2)'),
                                                "regensburg": gv.gvar('4.154(54)')      }
am_V = gv.gvar('1.41609(37)')
amp = gv.gvar('0.16381(54)')

## onemm - only at fine a 

plym_onemm = [gv.gvar('2.001(52)')*ainv['f-5']]     # 4-by-4 fit    
print(plym_onemm)
#plym_onemm_gev = [ plym_onemm[i]*(hbarc/a[i-1]) for i in range(len(plym_onemm)) ]
onemm_mass = {"hadspec"   : gv.gvar('4.411(30)'),
              "brambilla" : gv.gvar('4.410(15)'),
              "chen"      : gv.gvar('4.33(2)')    }

##Experimental states
onemm_exp = {"psi4230" : gv.gvar('4.220(15)'),
             "psi4360" : gv.gvar('4.368(13)'),
             "psi4415" : gv.gvar('4.421(4)')     }

