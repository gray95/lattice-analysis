import matplotlib.pyplot as plt
import numpy as np
import gvar as gv
from scipy.optimize import curve_fit
import sys
sys.path.append('..')
import math
from results import vc_up_diff, c_up_diff, f_up_diff
from results import vc_down_diff, c_down_diff, f_down_diff
from commonweal import w0overa, w0, hbarc
from gvar import log, sqrt

ydata  = np.concatenate((vc_up_diff[:-1], c_up_diff[:-1], f_up_diff[:-1]))
mq      = [ 3 , 5 , 7 ] 
A = { x:(w0/(hbarc*w0overa[x])).mean for x in ['vc','c','f'] }
xdata = [ (A[astr]**2,m) for astr in ['vc','c','f'] for m in mq ] 
print(xdata)

y     = [ t.mean for t in ydata   ] 
y_err = [ t.sdev for t in ydata   ] 

#plt.errorbar(m_l,y,y_err, fmt="ro")

#sys.exit(0)

def linear(x, c, m):
     return m*x + c

def multilinear(x, c, m1, m2):
    out = []
    for x in x:
        out.append(m1*x[0] + m2*x[1] + c)
    return out

def err_multilinear(x, c,c_err, m1,m1_err, m2,m2_err) :
    out = []
    for x in x:
        out.append(math.sqrt( x[0]**2*m1_err**2 + x[1]**2*m2_err**2 + c_err**2 ))
    return out

def err_linear(x,c,c_err, m,m_err) :
    return math.sqrt( x**2*m_err**2 + c_err**2 )

#popt, pcov = curve_fit(linear, xdata=m_l, ydata=y, sigma=y_err, p0 = [1, 1] )
popt, pcov = curve_fit(multilinear, xdata=xdata, ydata=y, sigma=y_err, p0 = [1, 1, 1] )
perr = np.sqrt(np.diag(pcov))

mm = np.linspace(0.5, 8, 50)    # scan over quark masses
a2 = np.linspace(0, A['vc']**2, 50)# scan over lattice spacings in GeV^-1
xx = [ (a,m) for a in a2 for m in mm ]
f_model = multilinear(xx,popt[0],popt[1],popt[2]) 

f_err = err_multilinear(xx,popt[0],perr[0],popt[1],perr[1],popt[2],perr[2])  

f_plus = f_model + f_err
f_minus = f_model + -1*f_err

print("Fit parameters c = %e" % popt[0] , "+/- %e" % perr[0] )
print("Fit parameters m1 = %e" % popt[1] , "+/- %e" % perr[1] )
print("Fit parameters m2 = %e" % popt[2] , "+/- %e" % perr[2] )

qmass = 1.0 
acont = 0.0
phys = [(0.0,1.0)]
f_ans     = multilinear(phys,popt[0],popt[1],popt[2])
f_ans_err = err_multilinear(phys,popt[0],perr[0],popt[1],perr[1],popt[2],perr[2])

print("Value at phys point is %s +- %s"%(f_ans, f_ans_err))

#plt.plot(xx, f_model, "r--")
#plt.xlim(0,9)

#plt.fill_between(xx, f_plus, f_minus, 
#                 color='gray', alpha=0.2)


#plt.xlabel(r'$\frac{m_q}{m_l}$', fontsize=12)
#plt.ylabel(r'$\delta a_{\mu}^{(\mathrm{light})}$' )

#plt.savefig("scipy_coarse.png")
#plt.show()
