import gvar as gv
import numpy as np
import sys
###[ 3ml , 5ml, 7ml ]  [amu, delta amu, ratio]

pickle23 = {"encoding":"latin1"}    # for compatibility btwn python2/3

vc_phys15= gv.load('../store/vcoarse/mphys15.p', **pickle23)
vc_phys= gv.load('../store/vcoarse/mphys.p', **pickle23)
vc_3ml = gv.load('../store/vcoarse/3ml.p', **pickle23)
vc_5ml = gv.load('../store/vcoarse/5ml.p', **pickle23)
vc_7ml = gv.load('../store/vcoarse/7ml.p', **pickle23)
vc_ms = gv.load('../store/vcoarse/ms_vcoarse.p', **pickle23)

c_3ml = gv.load('../store/coarse/3ml.p', **pickle23)
c_5ml = gv.load('../store/coarse/5ml.p', **pickle23)
c_7ml = gv.load('../store/coarse/7ml.p', **pickle23)
c_ms = gv.load('../store/coarse/ms_rho_coarse.p', **pickle23)

f_3ml = gv.load('../store/fine/3ml.p', **pickle23)
f_5ml = gv.load('../store/fine/5ml.p', **pickle23)
f_7ml = gv.load('../store/fine/7ml.p', **pickle23)
f_ms = gv.load('../store/fine/ms.p', **pickle23)

#############################################################
################ VERY COARSE - 0.15fm #######################

vc_amu = { "mphys":{1.5:[], 2:[]},"3ml":{1.5:[], 2:[]},"5ml":{1.5:[], 2:[]},"7ml":{1.5:[], 2:[]},"ms":{2.5:[]} }

vc_up_noqed    = vc_phys["up-nocharge"]
vc_up_phys_diff= vc_phys["up-charge"] - vc_phys["up-nocharge"]
vc_down_noqed  = vc_phys["down-nocharge"]
vc_down_phys_diff   = vc_phys["down-charge"] - vc_phys["down-nocharge"]
vc_amu_phys    = vc_phys["up-nocharge"]+vc_phys["down-nocharge"]
vc_amuqed_phys = vc_phys["up-charge"]+vc_phys["down-charge"]

vc_rt_phys     = vc_amuqed_phys/vc_amu_phys
vc_diff_phys   = vc_amuqed_phys - vc_amu_phys

vc_amu["mphys"][2].extend([vc_amu_phys, vc_diff_phys, vc_rt_phys])
vc_amu["3ml"][2].extend( vc_3ml["mq"] )
vc_amu["5ml"][2].extend( vc_5ml["mq"] )
vc_amu["7ml"][2].extend( vc_7ml["mq"] )
vc_amu["ms"][2.5].extend([vc_ms["amus"], vc_ms["diff"], vc_ms["ratio"]])

vc_amu_up = [x['up-nocharge'] for x in [vc_3ml, vc_5ml, vc_7ml] ]
vc_amu_down = [x['down-nocharge'] for x in [vc_3ml, vc_5ml, vc_7ml] ]
vc_up_diff = [x['up-charge']-x['up-nocharge'] for x in [vc_3ml, vc_5ml, vc_7ml] ]
vc_down_diff = [x['down-charge']-x['down-nocharge'] for x in [vc_3ml, vc_5ml, vc_7ml] ]
vc_diff = np.add(vc_up_diff, vc_down_diff)
#----------------------------------------#
vc_phi_m       = gv.gvar('1.0365(56)')
vc_phi_m_rt    = gv.gvar('1.000236(21)')
vc_phi_m_diff  = gv.gvar('0.245(22)') 	# MeV
vc_phi_f       = gv.gvar('241.0(1.8)')
vc_phi_f_rt    = gv.gvar('0.99984(12)')
vc_phi_f_diff  = gv.gvar('-0.038(28)')
## Q=2e/3 - unphysical
#vc_amu_rt_s    = gv.gvar('')
#vc_amudiff_s   = gv.gvar('')
#----------------------------------------#
#vc_phi_mass   = gv.gvar('')
#vc_delta_phi  = gv.gvar('')
#vc_fphi       = gv.gvar('')
#vc_fphi_diff    = gv.gvar('')

#############################################################
################ COARSE - 0.12fm #######################

c_up_diff = [x['up-charge']-x['up-nocharge'] for x in [c_3ml, c_5ml, c_7ml] ]
c_down_diff = [x['down-charge']-x['down-nocharge'] for x in [c_3ml, c_5ml, c_7ml] ]

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
f_up_diff = [x['up-charge']-x['up-nocharge'] for x in [f_3ml, f_5ml, f_7ml] ]
f_down_diff = [x['down-charge']-x['down-nocharge'] for x in [f_3ml, f_5ml, f_7ml] ]

#---------------------------------------#
f_phi_m       = gv.gvar('1.0268(57)')
f_phi_m_rt    = gv.gvar('1.0003(13)')
f_phi_m_diff  = gv.gvar('0.3(1.3)') 	# MeV
f_phi_f       = gv.gvar('238.0(2.2)')
f_phi_f_rt    = gv.gvar('1.0001(84)')
f_phi_f_diff  = gv.gvar('0.02(1.44)')
###########################################################


