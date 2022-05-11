import gvar as gv

## add in paths to corr files - need to be in gpl format

base = '/home/gray/Desktop/lattice-analysis/data/qqed'

vt_fname = {'vcphys':'vcoarse/mud_vt_vcoarse.gpl', 'vc':'vcoarse/vt_vcoarse.gpl', 'c':'coarse/vt_coarse.gpl', 'f':'fine/vt_fine.gpl'}
ps_fname = {'vc':'vcoarse/ps_vcoarse.gpl', 'c':'coarse/ps_coarse.gpl', 'f':'fine/ps_fine.gpl'}


hbarc = 0.197326968

w0 = gv.gvar('0.1715(9)')	# fm
w0_bmw = gv.gvar('0.17236(70)') # including IB

## Ensemble info and Renormalisation factors (for local currents)
w0overa = { "vc": gv.gvar('1.13215(35)'), 
            "c": gv.gvar('1.4149(6)'),
	    "f": gv.gvar('1.9518(7)'),
            "cFV": gv.gvar('1.4029(9)') }    #PRD 101,034512 (2020) Table I



## all numbers at scale of 2 GeV
ZV = {"vc": gv.gvar('0.95932(18)'), 
      "c": gv.gvar('0.97255(22)'),
      "f": gv.gvar('0.98445(11)'),
      "cFV": gv.gvar('0.97255(22)') } 

ZVqedu = {"vc": gv.gvar('0.999544(14)'),
	 "c": gv.gvar('0.999631(24)'),
	 "f": gv.gvar('0.999756(32)'),
         "cFV": gv.gvar('0.999631(24)')     } # Table X [2005.01845] Q=2/3  

ZVqedd = {}

for ensemble in ZVqedu:
    ZVqedd[ensemble] = 1 - (1-ZVqedu[ensemble])/4 

# EM quark renorm factor - following MILC - implies Dashen
delta_u = {'vc':gv.gvar('0.002049(77)'), 'c':gv.gvar('0.00208(12)'), 'f':gv.gvar('0.00210(50)')} 
delta_d = {'vc':gv.gvar('0.000503(16)'), 'c':gv.gvar('0.00053(10)'), 'f':gv.gvar('0.00055(50)')}

# HPQCD results [For summary: Phys Rev D 101, 034512 (2020) - Table VI]
hpqcd_amus = gv.gvar('53.41(59)e-10') 		# from 2014. isospin symm, no qed

# BMW
bmw_amus = gv.gvar('53.393(89)e-10')	# isospin symm, no qed, only stat error

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


