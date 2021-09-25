import gvar as gv
import sys
## Ensemble info and Renormalisation factors (for local currents)
hbarc = 0.197326968

w0 = gv.gvar('0.1715(9)')	# fm
a_vcoarse = w0/gv.gvar('1.13215(35)')
a_coarse = w0/gv.gvar('1.4029(9)')
a_fine = w0/gv.gvar('1.9006(20)')

a = [a_vcoarse, a_coarse, a_fine]
plym_onemp = [gv.gvar('3.292(83)'), gv.gvar('2.775(82)'), gv.gvar('1.934(80)')]
plym_onemp_gev = [ plym_onemp[i]*(hbarc/a[i]) for i in range(len(a)) ] 
parity_partner = [gv.gvar('2.99(19)'),gv.gvar('2.40(12)'),gv.gvar('1.53(22)')]
pp_gev = [ parity_partner[i]*(hbarc/a[i]) for i in range(len(a)) ] 

onemp_a = {"hadspec": 0.12, "ma": 0.18, "regensburg": 0.1145}

onemp_mass = {"hadspec": gv.gvar('4.310(29)'), 
						"ma": gv.gvar('4.309(2)'),
						"regensburg": gv.gvar('4.154(54)')	}

ZV=gv.gvar('0.977214(35)')
am_V = gv.gvar('1.41609(37)')
amp = gv.gvar('0.16381(54)')

## onemm - only at fine a 

plym_onemm = [gv.gvar('2.000(80)')]		# placeholder
plym_onemm_gev = [ plym_onemm[i]*(hbarc/a[i-1]) for i in range(len(plym_onemm)) ]
onemm_mass = {"hadspec" : gv.gvar('4.411(30)'),
							"brambilla" : gv.gvar('4.410(15)')	}

print(plym_onemm_gev)


