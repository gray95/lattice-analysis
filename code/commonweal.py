import gvar as gv
import sys
import os
basedir = '../data/hybrid'

## Ensemble info and Renormalisation factors (for local currents)

vc_phys = 'l3248f211b580m002426m06730m8447'
c_phys  = 'l4864f211b600m001907m05252m6382'
f_phys  = 'l6496f211b630m0012m0363m432'

c_unphys10 = 'l3264f211b600m00507m0507m628a-coul-v5'  
f_unphys5  = 'l3296f211b630m0074m037m440-coul-v5'

ensembles = {'vc':vc_phys, 'c':c_phys, 'c-10':c_unphys10, 'f':f_phys, 'f-5':f_unphys5} 

hbarc = 0.197326968

w0 = gv.gvar('0.1715(9)')	# fm
## Ensemble info and Renormalisation factors (for local currents)
w0overa = { "vc": gv.gvar('1.13215(35)'),
            "c": gv.gvar('1.4149(6)'),
            "c-10":gv.gvar('1.4029(9)'),
            "f": gv.gvar('1.9518(7)'),
            "f-5":gv.gvar('1.9006(20)') }    #PRD 101,034512 (2020) Table I

ZV=gv.gvar('0.977214(35)')  # ?

e_str = ['vc', 'c', 'c-10', 'f', 'f-5']
a_fm = { e:w0/w0overa[e] for e in e_str}
ainv = { e:hbarc/a_fm[e] for e in a_fm }
## paths to 1-+ correlators
print(a_fm)
onemp = {'vc':'onemp_vc_abij_m863.gpl', 'c':'onemp_coarsephys_m643a-h.gpl', 'c-10':'16tav_onemp_coarse_m650.gpl', 'f':'onemp_finephys_bph_-m432.gpl', 'f-5':'onemp_fine_m450.gpl'}

onemp_paths = { x:os.path.join(basedir,ensembles[x],onemp[x]) for x in ensembles } 

onemm_paths = {'f-5':os.path.join(basedir, ensembles['f-5'], 'tav_onemmz_m0450.gpl')}

