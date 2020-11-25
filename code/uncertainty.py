# A short (hopefully) script to compute uncertainties on various lattice quantities

import math as m
import numpy as np


Zv = 0.977214				# vector renorm constant
OZv= 0.000035				# std of vec renorm const

Z = Zv*0.197				# 0.197 for conversion to physical units

a  = 0.0886					# lattice spacing in fm
Oa = 0.0009					# std of a 

# vector ground state mass in lattice units and std
m0 = [1.41731, 1.4128, 1.42427, 1.41956, 1.41670, 1.4163, 1.4171, 1.4102, 1.4157]
Om0= [0.00060, 0.0031, 0.00096, 0.00096, 0.00098, 0.0020, 0.0020, 0.0043, 0.0041]			

# vector ground state amplitude in lattice units and std
a0 = [0.1650, 0.1566, 0.1743, 0.1682, 0.1646, 0.1643, 0.1659, 0.1566, 0.1631]				
Oa0= [0.0013, 0.0062, 0.0013, 0.0014, 0.0022, 0.0030, 0.0029, 0.0074, 0.0072]

for i in range(len(a0)):
			

	mphys = m0[i] * (0.197/a) 

	var_mphys = mphys**2 * ((Om0[i]/m0[i])**2 + (Oa/a)**2)

	print(mphys, m.sqrt(var_mphys)) 

	fv = Z * (a0[i]/a) * m.sqrt((2/m0[i])) * 10**3		# in MeV 

	var_fv = fv**2 * ((Oa0[i]**2/a0[i]**2) + (Oa**2/a**2) + (Om0[i]**2/(4*m0[i]**2))) 

	print(fv, m.sqrt(var_fv))
