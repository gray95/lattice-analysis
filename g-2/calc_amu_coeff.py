
#  Compute the w factors from RBC/UKQCD
#  arXiv:1512.09054
#

# maximum timeslice
T = 30 
##m_mu_gev = 0.106  
m_mu_gev = 0.1056583
alpha = 1/137.0 
##a_fm = 0.15
##ainv = 0.1973 / a_fm 
##ainv = 1.3026 

a_fm = 0.09
ainv = 0.1973 / a_fm 


print ("Mass muon = " , m_mu_gev , "  GeV")
##print ("Lattice spacing a = " , a_fm , "fm" )
print ("1/a = " , ainv , " GeV^-1")

m_mu = m_mu_gev  / ainv 

##m_mu = 1000 * m_mu_gev
##m_mu =  m_mu_gev

import math 
import pickle

#  https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html
from scipy.integrate import quad
import numpy as np

#
#  Equation 5 in the ukqcd/rbc paper  arXiv:1512.09054
#
def kernel(t,qsq) :
   q = math.sqrt(qsq)
   ans  = ((math.cos(q*t) - 1.0) / (qsq) + t*t/2 )
   return ans

#
#  f(K^2) from equation 3 in Blum's paper  hep-lat/0212018
#

def f(Ksq):
    Z = -(Ksq - math.sqrt(Ksq*Ksq + 4 * m_mu**2 * Ksq))/ ( 2.0 * m_mu**2 * Ksq )
    top = m_mu**2 * Ksq * Z**3 * ( 1 - Ksq * Z) 
    bot = 1 + m_mu**2 * Ksq * Z**2

    ans = top / bot
    return ans
 
#
#  integrand is q**2
#
def integrand(qsq,tt): 
    ans = kernel(tt,qsq) * f(qsq)
    return ans 



print ("Creation of w(t) factors for  a_mu")

print ("Maximum time = " , T)
print ("Mass of muon = " ,  m_mu , "  Gev")

print ("Test Kernel " ,  kernel(2,0.9) )
print ("Test Blum f " , f(0.8) )

t = 1 

print ("Test t = " , t , " Integrand =  " ,  integrand(0.8, 1) )



##  args=(a,b)
tmax = 30

wcoeff = np.zeros(tmax) 
tt     = np.zeros(tmax) 

for t in range(0,tmax) :
  result = quad(integrand, 0, np.inf, args=(t))
  rr = 4.0 * alpha**2 * result[0] 

  rr *= 2.0
##  rr  /= 9
##  print (t, result , rr)
  print (t, rr)
  wcoeff[t] = rr
  tt[t] = t 

fname="w_coeff"

print ("Saving coefficients to " , fname)


# Save to file
#np.save(fname, wcoeff)

##
##  plotting coefficients 
##

import matplotlib.pyplot as plt

plt.title("Weighted disconnected vector correlator ")
#plt.legend()

plt.xlabel('t')
plt.ylabel(r'$w(t)$')

plt.plot(tt, wcoeff)

#plt.savefig("wt.pdf")

##plt.xlim(0,20)

plt.show()
