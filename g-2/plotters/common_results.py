import gvar as gv

## OUR RESULTS - coarse
# qcd+qed values
c_amu = [gv.gvar('5.31(13)'), gv.gvar('4.872(77)'), gv.gvar('4.549(58)')]
# rt of qcdqed/qcd
c_rt = [gv.gvar('0.9980(12)'), gv.gvar('0.99845(54)'), gv.gvar('0.99862(31)')]
# diff - qcdqed-qcd/qcd
c_diff = [gv.gvar('-0.0020(12)'), gv.gvar('-0.00155(54)'), gv.gvar('-0.00138(31)')]
## VERY COARSE
vc_amu = [gv.gvar('5.09(10)'), gv.gvar('4.727(68)'), gv.gvar('4.416(54)')]
vc_rt  = [gv.gvar('0.9997(12)'), gv.gvar('0.99895(51)'), gv.gvar('0.99880(27)')]
vc_diff= [gv.gvar('-0.0003(12)'), gv.gvar('-0.00105(51)'), gv.gvar('-0.00120(27)')]
vc_amu_phys = gv.gvar('5.49(18)')		# combined a and b stream
vc_rt_phys  = gv.gvar('1.0057(65)')
vc_diff_phys= gv.gvar('0.0057(65)')


# unless specified ONLY qed isospin breaking included.

# BMW numbers - continuum values - only statistical errors included
bmw_amuqcd = gv.gvar('633.7(2.1)')	
bmw_amuqcdqed = bmw_amuqcd + gv.gvar('-1.23(40)')# + gv.gvar('6.60(63)')
bmw_amuqcdqedsib = bmw_amuqcd + gv.gvar('-1.23(40)') + gv.gvar('6.60(63)')
bmw_rt = bmw_amuqcdqed/bmw_amuqcd
bmw_rtsib = bmw_amuqcdqedsib/bmw_amuqcd	# including strong isospin breaking
bmw_diff = (bmw_amuqcdqed-bmw_amuqcd)/bmw_amuqcd
bmw_diffsib = (bmw_amuqcdqedsib-bmw_amuqcd)/bmw_amuqcd

# Giusti et al [1909.01962]
giu_diff = gv.gvar('1.1(1.0)')/gv.gvar('619.0(17.8)')
# note: they got the denominator from somewhere else

# HPQCD [For summary: Phys Rev D 101, 034512 (2020) - Table VI]
hpqcd_amu_ud = gv.gvar('637.8(8.8)e-10') # isospin symm, no qed



## STRANGE MASS
#----------------------------------------------
# Our results [VERY Coarse]- incl qQED
vc_amu_s = gv.gvar('53.70(58)e-10')
vc_amudiff_s = gv.gvar('-6.01(17)e-12')
vc_amudiffrt_s = gv.gvar('-0.001118(28)')
# [Coarse]
c_amu_s = gv.gvar('5.395(54)e-09')
c_amudiff_s = gv.gvar('-6.43(17)e-12')
c_amudiffrt_s = gv.gvar('-0.001190(29)')

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
