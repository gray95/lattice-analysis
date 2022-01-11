import gvar as gv


hbarc = 0.197326968

w0_bmw = gv.gvar('0.17236(70)') # including IB

## Ensemble info and Renormalisation factors (for local currents)
w0overa = { "vc": gv.gvar('1.13215(35)'), 
            "c": gv.gvar('1.41490(60)'),
	    "f": gv.gvar('1.95180(70)')	} #PRD 101,034512 (2020) Table I

w0 = gv.gvar('0.1715(9)')	# fm

ZV = {"vc": gv.gvar('0.9881(10)'), 
      "c": gv.gvar('0.99220(40)'),
      "f": gv.gvar('0.99400(50)')		} # at m_s

ZVqed = {"vc": gv.gvar('0.999544(14)'),
	 "c": gv.gvar('0.999631(24)'),
	 "f": gv.gvar('0.999756(32)')     } # at m_l, PRD 100,114513 (2019) Table IV 



