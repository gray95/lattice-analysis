from math import pi
from math import sqrt
import gvar as gv
from commonweal import hbarc, ZV, a_fine #, amp, am_V

am_V = gv.gvar('1.3659(49)')
amp  = gv.gvar('0.1895(54)') 	# vector operator INTO hybrid state

alpha = 1/134 		# qed coupling at charm mass

m_V = hbarc/a_fine * am_V	* 10**3	# MeV

f_V = ZV*hbarc*amp/a_fine * (2/am_V)**0.5 * 10**3 # MeV

vec_lep_width = (16*pi/27) * alpha**2 * f_V**2 * (1/m_V) * 10**3		# all in MeV

print("mass is ", m_V, "MeV")
print("vector decay constant is ", f_V, "MeV")
print("vector leptonic decay width is: ", vec_lep_width, "keV")




## inverse problem ( decay constant from leptonic width) any vector state should work
exp_vec_lep_width = gv.gvar('2.5(3)') * 10**(-6) 		# MeV

_m_V = gv.gvar('4230(100)')
_f_V = gv.sqrt( exp_vec_lep_width * (27/(16*pi)) * (_m_V/alpha**2) )

print ("exotic psi(4230) vector decay constant is ", _f_V, "MeV" )
