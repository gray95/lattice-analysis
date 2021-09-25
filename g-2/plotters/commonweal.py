import gvar as gv


hbarc = 0.197326968



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
				 "f": gv.gvar('0.999756(32)')		} # at m_l, PRD 100,114513 (2019) Table IV 


### LIGHT QUARK MASS ( << ms ) [ 3ml , 5ml, 7ml ]
## OUR RESULTS - coarse t*~2fm
# qcd+qed values
c_amu = [gv.gvar('5.188(98)e-08'), gv.gvar('4.826(62)e-08'), gv.gvar('4.529(50)e-08')]
c_rt = [gv.gvar('0.9987(31)'), gv.gvar('0.9985(14)'), gv.gvar('0.99862(79)')]
c_diff  = [gv.gvar('-7(16)e-11'), gv.gvar('-7.2(7.0)e-11'), gv.gvar('-6.3(3.6)e-11')]
c_diff_rt = [gv.gvar('-0.0013(31)'), gv.gvar('-0.0015(14)'), gv.gvar('-0.00138(79)')]
c_diff_phys_ext = gv.gvar('-6.4(4.5)e-11')
c_rt_phys_ext = gv.gvar('0.9985(33)')

## VERY COARSE t*~2fm
vc_amu = [gv.gvar('5.158(91)e-08'), gv.gvar('4.799(62)e-08'), gv.gvar('4.497(51)e-08')]
vc_rt  = [gv.gvar('0.9994(29)'), gv.gvar('0.9990(16)'), gv.gvar('0.99889(86)')]
vc_diff_rt= [gv.gvar('-0.0006(29)'), gv.gvar('-0.0010(16)'), gv.gvar('-0.00111(86)')]
vc_diff   = [gv.gvar('-3(15)e-11'), gv.gvar('-4.7(7.5)e-11'), gv.gvar('-5.0(3.9)e-11')]
vc_rt_phys_ext = gv.gvar('0.9995(34)')


vc_amu_phys = gv.gvar('5.47(22)e-08')		# combined b streams - none of dan's a stream
vc_rt_phys  = gv.gvar('1.0155(81)') # incl sib
vc_diff_phys= gv.gvar('0.0155(81)')
vc_diff_phys_ext = gv.gvar('-4.7(4.1)e-11')


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



## STRANGE MASS t* = 2.5fm
#----------------------------------------------
# Our results [VERY Coarse]- incl qQED
vc_amu_s       = gv.gvar('5.463(56)e-09')
vc_amu_rt_s    = gv.gvar('0.998877(29)')
vc_amudiff_s   = gv.gvar('-6.14(17)e-12')
vc_amudiffrt_s = gv.gvar('-0.001123(29)')
#----------------------------------------#
vc_phi_mass   = gv.gvar('1.0382(56)')
vc_delta_phi  = gv.gvar('0.256(19)')
vc_fphi       = gv.gvar('0.1580(11)')
vc_fphi_rt    = gv.gvar('0.999994(89)')
# [Coarse]
c_amu_s       = gv.gvar('5.386(54)e-09')
c_amu_rt_s    = gv.gvar('0.998989(48)')
c_amudiff_s   = gv.gvar('-5.45(27)e-12')
c_amudiffrt_s = gv.gvar('-0.001011(48)')
#---------------------------------------#
c_phi_mass    = gv.gvar('1.0332(57)')
c_delta_phi   = gv.gvar('0.26(64)') 	# MeV
c_fphi        = gv.gvar('0.1381(12)')
c_fphi_rt     = gv.gvar('1.0001(32)')
# [Fine]
f_amu_s       = gv.gvar('5.356(55)e-09')
f_amu_rt_s    = gv.gvar('1.0057(98)')
f_amudiff_s   = gv.gvar('3.0(5.2)e-11')
f_amudiffrt_s = gv.gvar('0.0057(98)')
#---------------------------------------#
f_phi_mass    = gv.gvar('1.034(11)')
f_delta_phi   = gv.gvar('-6.6(9.9)') 	# MeV
f_fphi        = gv.gvar('0.1190(36)')
f_fphi_rt     = gv.gvar('0.966(32)')
###########################################################

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
