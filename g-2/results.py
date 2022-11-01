import gvar as gv
import numpy as np
import sys

###                         [ 3ml ,5ml, 7ml, ms ] 

vc = gv.load('/home/gray/Desktop/lattice-analysis/g-2/out/vc_nofit.p')
c = gv.load('/home/gray/Desktop/lattice-analysis/g-2/out/c_nofit.p')
f = gv.load('/home/gray/Desktop/lattice-analysis/g-2/out/f_nofit.p')

#vc = gv.load('/home/gray/Desktop/lattice-analysis/g-2/fits/vt_vc_bin2_tmin2.p')
#c = gv.load('/home/gray/Desktop/lattice-analysis/g-2/fits/vt_c_bin2_tmin2.p')
#f = gv.load('/home/gray/Desktop/lattice-analysis/g-2/fits/vt_f_bin2_tmin2.p')


# compute diffs
vc_amud      = vc['down-nocharge']
vc_down_diff = [ (y2-y1) for (y2,y1) in zip(vc['down-charge'],vc['down-nocharge']) ]
vc_up_diff =   [ (y2-y1) for (y2,y1) in zip(vc['up-charge'],vc['up-nocharge']) ]

c_amud      = c['down-nocharge']
c_down_diff = [ (y2-y1) for (y2,y1) in zip(c['down-charge'],c['down-nocharge']) ]
c_up_diff =   [ (y2-y1) for (y2,y1) in zip(c['up-charge'],c['up-nocharge']) ]

f_amud      = f['down-nocharge']
f_down_diff = [ (y2-y1) for (y2,y1) in zip(f['down-charge'],f['down-nocharge']) ]
f_up_diff =   [ (y2-y1) for (y2,y1) in zip(f['up-charge'],f['up-nocharge']) ]

#############################################################

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

