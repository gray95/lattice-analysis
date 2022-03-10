import matplotlib.pyplot as plt
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
from scipy.optimize import curve_fit
import sys
import math

c_diff  = [gv.gvar('-1.8(3.0)e-11'), gv.gvar('-3.3(1.1)e-11'), gv.gvar('-4.06(62)e-11')]
m_l      = [ 3 , 5 , 7 ] 

if 0 :
    print(c_diff[0].mean)
    print(c_diff[0].sdev)

y     = [ t.mean for t in c_diff   ] 
y_err = [ t.sdev for t in c_diff   ] 

print(y)
print(y_err)
print(m_l)

plt.errorbar(m_l,y,y_err, fmt="ro")



##sys.exit(0)

def linear(x, c, m):
     return m*x + c

def err_linear(x,c,c_err, m,m_err) :
    return math.sqrt( x**2*m_err**2 + c_err**2 )

popt, pcov = curve_fit(linear, xdata=m_l, ydata=y, sigma=y_err, p0 = [1, 1] )

perr = np.sqrt(np.diag(pcov))

xx = np.linspace(0, 8, 50)
f_model = np.array([linear(f,popt[0],popt[1]) for f in xx])

f_err = np.array( [ err_linear(f,popt[0],perr[0],popt[1],perr[1]) for f in xx]     )  

f_plus = f_model + f_err
f_minus = f_model - f_err


print("Fit parameters c = %e" % popt[0] , "+/- %e" % perr[0] )
print("Fit parameters m = %e" % popt[1] , "+/- %e" % perr[1] )

qmass = 1.0 
f_ans     = linear(qmass,popt[0],popt[1])
f_ans_err = err_linear(qmass,popt[0],perr[0],popt[1],perr[1])

print("Value at ml = 1 = ", f_ans , "+/-" , f_ans_err)

plt.plot(xx, f_model, "r--")
plt.xlim(0,9)

plt.fill_between(xx, f_plus, f_minus, 
                 color='gray', alpha=0.2)


plt.xlabel(r'$\frac{m_q}{m_l}$', fontsize=12)
plt.ylabel(r'$\delta a_{\mu}^{(\mathrm{light})}$' )

#plt.savefig("scipy_coarse.png")
plt.show()
