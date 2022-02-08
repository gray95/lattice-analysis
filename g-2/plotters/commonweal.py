import gvar as gv

## add in paths to corr files - need to be in gpl format

f_ms  = "../data/qqed/fine/ms_fine.gpl"
c_ms  = ""
vc_ms = ""

hbarc = 0.197326968

w0_bmw = gv.gvar('0.17236(70)') # including IB

## Ensemble info and Renormalisation factors (for local currents)
w0overa = { "vc": gv.gvar('1.13215(35)'), 
            "c": gv.gvar('1.4149(6)'),
	    "f": gv.gvar('1.9518(7)'),
            "cFV": gv.gvar('1.4029(9)') }    #PRD 101,034512 (2020) Table I

w0 = gv.gvar('0.1715(9)')	# fm

## renormalisation ov loval vector current at ms?


ZV = {"vc": gv.gvar('0.95932(18)'), 
      "c": gv.gvar('0.97255(22)'),
      "f": gv.gvar('0.98445(11)'),
      "cFV": gv.gvar('0.97255(22)') } 

ZVqed = {"vc": gv.gvar('0.999544(14)'),
	 "c": gv.gvar('0.999631(24)'),
	 "f": gv.gvar('0.999756(32)'),
         "cFV": gv.gvar('0.999631(24)')     } # Table X [2005.01845] 


# HPQCD results [For summary: Phys Rev D 101, 034512 (2020) - Table VI]
hpqcd_amu_s = gv.gvar('53.41(59)e-10') 		# from 2014. isospin symm, no qed

# BMW
bmw_amu_s = gv.gvar('53.393(89)e-10')	# isospin symm, no qed, only stat error

# Giusti et al.
giu_amudiff_s = gv.gvar('-0.0053(33)e-10')
giu_amudiffrt_s = gv.gvar('-0.00010(7)') 

## CHARM MASS
#--------------------------------------------
# HPQCD 
hpqcd_amu_c = gv.gvar('14.40(40)e-10')
bmw_amu_c = gv.gvar('14.6(1)e-10')

# unless specified ONLY qed isospin breaking included.

# BMW numbers - continuum values - only statistical errors included
bmw_amuqcd = gv.gvar('633.7(2.1)')	
bmw_amuqcdqed = bmw_amuqcd + gv.gvar('-1.23(40)')# + gv.gvar('6.60(63)')
bmw_amuqcdqedsib = bmw_amuqcd + gv.gvar('-1.23(40)') + gv.gvar('6.60(63)')
bmw_rt = bmw_amuqcdqed/bmw_amuqcd
bmw_rtsib = bmw_amuqcdqedsib/bmw_amuqcd	# including strong isospin breaking
bmw_diff = gv.gvar('-1.23(40)')
bmw_diff_rt = (bmw_amuqcdqed-bmw_amuqcd)/bmw_amuqcd
bmw_diffsib = (bmw_amuqcdqedsib-bmw_amuqcd)/bmw_amuqcd

# Giusti et al [1909.01962]
giu_diff = gv.gvar('1.1(1.0)')/gv.gvar('619.0(17.8)')
# note: they got the denominator from somewhere else

# HPQCD [For summary: Phys Rev D 101, 034512 (2020) - Table VI]
hpqcd_amu_ud = gv.gvar('637.8(8.8)e-10') # isospin symm, no qed


