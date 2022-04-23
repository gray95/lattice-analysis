import gvar as gv
import numpy as np
import sys
###[ 3ml , 5ml, 7ml ]  [amu, delta amu, ratio]

pickle23 = {"encoding":"latin1"}    # for compatibility btwn python2/3

#vc = gv.load('../out/svd_test_vt_vc_bin1.p', **pickle23)
#c = gv.load('../out/svd_test_vt_c_bin1.p', **pickle23)
#f = gv.load('../out/svd_test_vt_f_bin1.p', **pickle23)
vc_ms_diff = gv.gvar('-2.20(15)e-12')
c_ms_diff = gv.gvar('-2.37(26)e-12')
f_ms_diff = gv.gvar('-2.54(48)e-12')

vc_up_diff = gv.gvar(['-2.7(3.1)e-11', '-3.3(1.2)e-11', '-3.40(59)e-11'])
c_up_diff = gv.gvar(['-3.7(9.8)e-11', '-3.6(4.3)e-11', '-3.8(2.3)e-11'])
f_up_diff = gv.gvar(['2(37)e-11', '-2(19)e-11', '-3(12)e-11'])

vc_down_diff = gv.gvar(['-1.6(4.3)e-12', '-2.0(2.2)e-12', '-2.1(1.1)e-12'])
c_down_diff = gv.gvar(['-3(24)e-12', '-2(11)e-12', '-2.4(5.8)e-12'])
f_down_diff = gv.gvar(['1(93)e-12', '-1(48)e-12', '-2(31)e-12'])
#############################################################
################ VERY COARSE - 0.15fm #######################

c_amus        = gv.gvar('5.173(52)e-09') # on physical ensemble
c_amus_rt     = gv.gvar('0.998990(48)')
c_amus_diff   = gv.gvar('-5.23(25)e-12')
#---------------------------------------#
c_phi_m       = gv.gvar('1.0336(55)')
c_phi_m_rt    = gv.gvar('1.000234(12)')
c_phi_m_diff  = gv.gvar('0.242(12)') 	# MeV
c_phi_f       = gv.gvar('240.2(1.6)')
c_phi_f_rt    = gv.gvar('0.999861(62)')
c_phi_f_diff  = gv.gvar('-0.033(15)')
# FV runs on L=24,40 unphysical - 1/10
c_24_amus        = gv.gvar('5.310(54)e-09')
c_24_amus_rt     = gv.gvar('0.999033(56)')
c_24_amus_diff   = gv.gvar('-5.14(31)e-12')
#---------------------------------------#
c_24_phi_m       = gv.gvar('1.0342(62) ')
c_24_phi_m_rt    = gv.gvar('1.0003(12)')
c_24_phi_m_diff  = gv.gvar('0.3(1.2)') # MeV
c_24_phi_f       = gv.gvar('242.6(2.7)') # MeV
c_24_phi_f_rt    = gv.gvar('1.0002(48)')
c_24_phi_f_diff  = gv.gvar('0.06(1.17)')
#--------------------------------------------#
c_40_amus        = gv.gvar('5.308(54)e-09')
c_40_amus_rt     = gv.gvar('0.999015(57)')
c_40_amus_diff   = gv.gvar('-5.23(31)e-12')
#---------------------------------------#
c_40_phi_m       = gv.gvar('1.0361(62)')
c_40_phi_m_rt    = gv.gvar('1.0002(21)')
c_40_phi_m_diff  = gv.gvar('0.3(2.2)') 	# MeV
c_40_phi_f       = gv.gvar('244.4(2.6)')
c_40_phi_f_rt    = gv.gvar('1.0001(79)')
c_40_phi_f_diff  = gv.gvar('0.02(1.93)')
## Q=2e/3 - unphysical - just the diffs
c_2Q_amus_diff   = gv.gvar('-9.48(27)e-12')
c_2Q_amus_rt     = gv.gvar('0.998170(49)')
c_2Q_phi_m_rt  = gv.gvar('1.000942(40)') 	# MeV
c_2Q_phi_f_rt  = gv.gvar('1.00103(20)')
#------------------------------------------#
c_24_2Q_amus_diff = gv.gvar('-8.57(49)e-12')
c_24_2Q_amus_rt   = gv.gvar('0.998387(90)')
c_24_2Q_phi_m_rt  = gv.gvar('1.0009(13)') 	# MeV
c_24_2Q_phi_f_rt  = gv.gvar('1.0010(53)')
#---------------------------------------------#
c_40_2Q_amus_diff = gv.gvar('-9.20(50)e-12')
c_40_2Q_amus_rt   = gv.gvar('0.998269(92)')
c_40_2Q_phi_m_rt  = gv.gvar('1.0010(21)') 	# MeV
c_40_2Q_phi_f_rt  = gv.gvar('1.0014(79)')

# order for FV lists is [24, 40, 48]
coarseFV_amus = [ c_24_amus, c_40_amus, c_amus ]
coarseFV_amus_diff = [ c_24_amus_diff, c_40_amus_diff, c_amus_diff ]
coarseFV_amus_rt = [ c_24_amus_rt, c_40_amus_rt, c_amus_rt ]  
coarseFV_phi_m_diff = [ c_24_phi_m_diff, c_40_phi_m_diff, c_phi_m_diff ]
coarseFV_phi_f_diff = [ c_24_phi_f_diff, c_40_phi_f_diff, c_phi_f_diff ]
coarseFV_phi_m_rt = [ c_24_phi_m_rt, c_40_phi_m_rt, c_phi_m_rt ]
coarseFV_phi_f_rt = [ c_24_phi_f_rt, c_40_phi_f_rt, c_phi_f_rt ]
coarseFV2Q_amus_diff = [ c_24_2Q_amus_diff, c_40_2Q_amus_diff, c_2Q_amus_diff ]
coarseFV2Q_amus_rt = [ c_24_2Q_amus_rt, c_40_2Q_amus_rt, c_2Q_amus_rt ]  
coarseFV2Q_phi_m_rt = [ c_24_2Q_phi_m_rt, c_40_2Q_phi_m_rt, c_2Q_phi_m_rt ]
coarseFV2Q_phi_f_rt = [ c_24_2Q_phi_f_rt, c_40_2Q_phi_f_rt, c_2Q_phi_f_rt ]

# [Fine]
#---------------------------------------#
f_phi_m       = gv.gvar('1.0268(57)')
f_phi_m_rt    = gv.gvar('1.0003(13)')
f_phi_m_diff  = gv.gvar('0.3(1.3)') 	# MeV
f_phi_f       = gv.gvar('238.0(2.2)')
f_phi_f_rt    = gv.gvar('1.0001(84)')
f_phi_f_diff  = gv.gvar('0.02(1.44)')
###########################################################


