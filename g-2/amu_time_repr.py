import sys
import os
import bz2 
import dill as pickle
#  https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html
from scipy.integrate import quad
import numpy as np
import gvar as gv
from gvar import log, sqrt, cos
from commonweal import w0, w0overa, ZV, ZVqedd, ZVqedu, hbarc
from commonweal import base, vt_fname, m_mu_gev, alpha, vtfitobj
import corrfitter as cf
import blinding

masses = ['3ml', '5ml', '7ml', 'ms']

## ensemble
a_str = 'vc'

a_fm = w0/w0overa[a_str]
ainv = hbarc / a_fm 
m_mu = m_mu_gev  / ainv 

print ("Mass muon = " , m_mu_gev , "  GeV")
print ("Lattice spacing a = " , a_fm , "fm" )
print ("1/a = " , ainv , " GeV^-1")

ZV = ZV[a_str]
ZVqedd = ZVqedd[a_str]*ZV
ZVqedu = ZVqedu[a_str]*ZV

B = 2
f = os.path.join(base, vt_fname[a_str])
corr = gv.dataset.Dataset(f)

##blinding
bF = "g2-conn-blindfactor.bin"
print("load blinding factor from %s"%bF)
blindF = blinding.read_blind(bF)

print("file %s \nstream %s\nbinsize = %d"%(f,a_str,B))
dset = gv.BufferDict()

# make sure corrs are correctly normalised
for tag in corr.keys():
    corr[tag] = 3*blindF*np.array(corr[tag])

corrbin = gv.dataset.bin_data(corr, binsize=B)
print("%d cfgs"%corrbin.samplesize)

data = gv.dataset.avg_data(corrbin)

keys = dict.fromkeys(masses)
for m in keys:
    keys[m] = {}
    K = [k for k in data.keys() if m in k]   
    keys[m]['n'] = K[0]
    keys[m]['d'] = K[1]
    keys[m]['u'] = K[2]

T = data[keys['3ml']['n']].shape[0]
print("lattice time extent: %d"%T)

#
#  Equation 5 in the ukqcd/rbc paper  arXiv:1512.09054
#
def kernel(t,qsq) :
   q = sqrt(qsq)
   ans = ( ((cos(q*t) - 1.0) / qsq) + t*t/2 ) 
   return ans

#
#  f(K^2) from equation 3 in Blum's paper  hep-lat/0212018
#

def f(Ksq):
#    Ksq = Ksq/ainv**2
    Z = -(Ksq - sqrt(Ksq*Ksq + 4 * m_mu**2 * Ksq))/ ( 2.0 * m_mu**2 * Ksq )
    top = m_mu**2 * Ksq * Z**3 * ( 1 - Ksq * Z) 
    bot = 1 + m_mu**2 * Ksq * Z**2
    ans = top / bot
    return ans
 
#
#  integrand is q**2
#
def integrand(qsq,t): 
    ans = 4.0 * alpha**2 * kernel(t,qsq) * f(qsq)
    return ans.mean

def amu_from_time_repr(G, tmin, tmax, T, Z=1., Q=1., periodic=None):
    amu = gv.BufferDict()
    if periodic:
        C = [ G[t]+G[-t] for t in range(T//2) ]
    else:
        C = [ 2*G[t] for t in range(T) ]

    for t in range(tmin,tmax):
        wt = quad(integrand, 0, np.inf, args=(t), epsabs=0)
        wt = gv.gvar(wt[0],wt[1])
        amut = wt * Z**2 * Q**2 * C[t] 
        amu[t] = amut
    return amu

## compute amu over different masses and charges
amudown = dict.fromkeys(masses) 
amuup = dict.fromkeys(masses) 
amudowndiff = dict.fromkeys(masses) 
amuupdiff = dict.fromkeys(masses) 

fit= pickle.load(bz2.BZ2File(vtfitobj[a_str], 'rb'))
tmin=0
tstar = T//2 - 1 #int(2/a_fm.mean)
print("t* is %d"%tstar)
tcut = T
for m in masses:
    print("mass = %s"%m)

    tags = ['no-charge', 'down-charge', 'up-charge']

    modeluncharged = blindF*cf.Corr2(datatag=tags[0],tdata=range(T),a=(m+':a:n',m+':ao:n'),b=(m+':a:n',m+':ao:n'),dE=(m+':dE:n',m+':dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargeddown = blindF*cf.Corr2(datatag=tags[1],tdata=range(T),a=(m+':a:d',m+':ao:d'),b=(m+':a:d',m+':ao:d'),dE=(m+':dE:d',m+':dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
    modelchargedup = blindF*cf.Corr2(datatag=tags[2],tdata=range(T),a=(m+':a:u',m+':ao:u'),b=(m+':a:u',m+':ao:u'),dE=(m+':dE:u',m+':dEo:u'),s=(1.,-1.)).fitfcn(fit.p)
    
    uncharged = np.append( data[keys[m]['n']][:tstar], modeluncharged[tstar:] ) 
    chargeddown = np.append( data[keys[m]['d']][:tstar], modelchargeddown[tstar:] ) 
    chargedup = np.append( data[keys[m]['u']][:tstar], modelchargedup[tstar:] ) 

    unchargedamud = amu_from_time_repr( uncharged, tmin=tmin, tmax=tcut, T=T, Z=ZV, Q=1./3, periodic=False )
    unchargedamuu = amu_from_time_repr( uncharged, tmin=tmin, tmax=tcut, T=T, Z=ZV, Q=2./3 , periodic=False)
    chargedamud   = amu_from_time_repr( chargeddown, tmin=tmin, tmax=tcut, T=T, Z=ZVqedd, Q=1./3 , periodic=False)
    chargedamuu   = amu_from_time_repr( chargedup, tmin=tmin, tmax=tcut, T=T, Z=ZVqedu, Q=2./3 , periodic=False)

    amud = amuu = qedamud = qedamuu = gv.gvar('0(0)')
    amudown[m] = gv.BufferDict()
    amuup[m] = gv.BufferDict()
    amudowndiff[m] = gv.BufferDict()
    amuupdiff[m] = gv.BufferDict()
    for t in range(tmin,T):
        amud += unchargedamud[t] 
        amuu += unchargedamuu[t] 
        qedamuu += chargedamuu[t] 
        qedamud += chargedamud[t] 
        downdiff = qedamud - amud
        updiff =   qedamuu - amuu

        amudown[m][t] = amud
        amuup[m][t]  = amuu
        amudowndiff[m][t] = downdiff
        amuupdiff[m][t] =   updiff
    print("total = %s"%(amud+amuu))
print(amudowndiff['ms'])


############################################################################
#   print("Switching to MODEL\n")

#   modeluncharged    = gv.BufferDict()
#   modelchargeddown  = gv.BufferDict()
#   modelchargedup    = gv.BufferDict()
#   model_amudown = dict.fromkeys(masses) 
#   model_amuup = dict.fromkeys(masses) 
#   model_amudowndiff = dict.fromkeys(masses) 
#   model_amuupdiff = dict.fromkeys(masses) 

#   tmin = 2
#   tmax = T

#   for mq in masses:
#       m=mq+':'
#       tags = ['no-charge', 'down-charge', 'up-charge']

#       modeluncharged[mq] = blindF*cf.Corr2(datatag=tags[0],tdata=range(T),a=(m+'a:n',m+'ao:n'),b=(m+'a:n',m+'ao:n'),dE=(m+'dE:n',m+'dEo:n'),s=(1.,-1.)).fitfcn(fit.p)
#       modelchargeddown[mq] = blindF*cf.Corr2(datatag=tags[1],tdata=range(T),a=(m+'a:d',m+'ao:d'),b=(m+'a:d',m+'ao:d'),dE=(m+'dE:d',m+'dEo:d'),s=(1.,-1.)).fitfcn(fit.p)
#       modelchargedup[mq] = blindF*cf.Corr2(datatag=tags[2],tdata=range(T),a=(m+'a:u',m+'ao:u'),b=(m+'a:u',m+'ao:u'),dE=(m+'dE:u',m+'dEo:u'),s=(1.,-1.)).fitfcn(fit.p)

#       model_unchargedamud = amu_from_time_repr( modeluncharged[mq], tmin=tmin, tmax=T, T=T, Z=ZV, Q=1./3 )
#       model_unchargedamuu = amu_from_time_repr( modeluncharged[mq], tmin=tmin, tmax=T, T=T, Z=ZV, Q=2./3 )
#       model_chargedamud   = amu_from_time_repr( modelchargeddown[mq], tmin=tmin, tmax=T, T=T, Z=ZVqedd, Q=1./3 )
#       model_chargedamuu   = amu_from_time_repr( modelchargedup[mq], tmin=tmin, tmax=T, T=T, Z=ZVqedu, Q=2./3 )

#       amud = amuu = qedamud = qedamuu = gv.gvar('0(0)')
#       model_amudown[mq] = gv.BufferDict()
#       model_amuup[mq] = gv.BufferDict()
#       model_amudowndiff[mq] = gv.BufferDict()
#       model_amuupdiff[mq] = gv.BufferDict()
#       for t in range(tmin,tmax):
#           amud += model_unchargedamud[t] 
#           amuu += model_unchargedamuu[t] 
#           qedamuu += model_chargedamuu[t] 
#           qedamud += model_chargedamud[t] 
#           downdiff = qedamud - amud
#           updiff =   qedamuu - amuu

#           model_amudown[mq][t] = amud
#           model_amuup[mq][t]  = amuu
#           model_amudowndiff[mq][t] = downdiff
#           model_amuupdiff[mq][t] =   updiff

## PLOTTING
import matplotlib.pyplot as plt

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':16})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)

fig, (ax1, ax2) = plt.subplots(2,1,sharex=True, figsize=(6,6))

Y =  amudowndiff['ms']
YY = amudowndiff['7ml']
U =  amuupdiff['ms']
UU = amuupdiff['7ml']

data_t = [t for t in Y if t<=tstar]
model_t = [t for t in Y if t>tstar]
data_t_fm = [ t*a_fm for t in data_t] 
model_t_fm = [ t*a_fm for t in model_t] 

data_y =  [Y[t] for t in data_t]
model_y = [Y[t] for t in model_t]
data_yy =  [YY[t] for t in data_t]
model_yy = [YY[t] for t in model_t]

ax1.errorbar( [t.mean for t in data_t_fm], [y.mean for y in data_y], xerr=[t.sdev for t in data_t_fm], yerr=[y.sdev for y in data_y], fmt='h', alpha=0.6, color='red', label='DATA' )
ax1.errorbar( [t.mean for t in model_t_fm], [y.mean for y in model_y], xerr=[t.sdev for t in model_t_fm], yerr=[y.sdev for y in model_y], fmt='h', alpha=0.6, color='black', label='MODEL' )
ax1.errorbar( [t.mean for t in data_t_fm], [y.mean for y in data_yy], xerr=[t.sdev for t in data_t_fm], yerr=[y.sdev for y in data_yy], fmt='h', alpha=0.6, color='red')
ax1.errorbar( [t.mean for t in model_t_fm], [y.mean for y in model_yy], xerr=[t.sdev for t in model_t_fm], yerr=[y.sdev for y in model_yy], fmt='h', alpha=0.6, color='black')

data_u =  [U[t] for t in data_t]
model_u = [U[t] for t in model_t]
data_uu =  [UU[t] for t in data_t]
model_uu = [UU[t] for t in model_t]

ax2.errorbar( [t.mean for t in data_t_fm], [y.mean for y in data_u], xerr=[t.sdev for t in data_t_fm], yerr=[y.sdev for y in data_u], fmt='h', alpha=0.6, color='red', label='DATA' )
ax2.errorbar( [t.mean for t in model_t_fm], [y.mean for y in model_u], xerr=[t.sdev for t in model_t_fm], yerr=[y.sdev for y in model_u], fmt='h', alpha=0.6, color='black', label='MODEL' )
ax2.errorbar( [t.mean for t in data_t_fm], [y.mean for y in data_uu], xerr=[t.sdev for t in data_t_fm], yerr=[y.sdev for y in data_uu], fmt='h', alpha=0.6, color='red')
ax2.errorbar( [t.mean for t in model_t_fm], [y.mean for y in model_uu], xerr=[t.sdev for t in model_t_fm], yerr=[y.sdev for y in model_uu], fmt='h', alpha=0.6, color='black')


ax1.set_xlim(left=1,right=4)
ax2.set_xlabel('t [fm]')
ax1.set_ylabel(r'$\Delta a_\mu$', rotation=0, labelpad=20, fontsize=20)
ax2.set_ylabel(r'$\Delta a_\mu$', rotation=0, labelpad=20, fontsize=20)
ax1.axvline(x=a_fm.mean*(T//2-1), color='black', linestyle='--', lw=1)
ax2.axvline(x=a_fm.mean*(T//2-1), color='black', linestyle='--', lw=1)
handles,labels = fig.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
ax1.legend(handles=handles,labels=labels,frameon=False,loc='upper center', fontsize=12)

plt.savefig('./figures/tstar_nofit_'+a_str+'.png', dpi=500, bbox_inches="tight")
plt.show()


