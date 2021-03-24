import matplotlib
matplotlib.use('Agg')
import numpy as np
import gvar as gv
from gvar import log,exp,evalcov
import lsqfit as lsq
from corrfitter import Corr2,Corr3,CorrFitter
import g2tools as g2
import matplotlib.pyplot as plt
import copy
import sys

# keep data files in one place
data_dir="../data/qed/vcoarse/dan_a"
pfile = "last_fit_tmp.bin" # last fit


def build_prior(nexp,dErho=gv.gvar('0.476(100)')):
    """ build prior """
    prior = gv.BufferDict()

    prior.add('a:rho',[gv.gvar(0,0.1),gv.gvar(0,0.1),gv.gvar(0,0.1),gv.gvar(0,0.1)])
    # fix it so that the rho terms are opposite sign to the omega
    prior.add('b:rho',[-prior['a:rho'][0],prior['a:rho'][1],prior['a:rho'][2],-prior['a:rho'][3]])
    prior.add('log(dE:rho)',[log(gv.gvar(0.5,0.4)) for i in range(nexp)])
    prior['log(dE:rho)'][0] = log(dErho) # rho
    prior['log(dE:rho)'][1] = log(gv.gvar('0.09(18)')) # omega
    prior['log(dE:rho)'][2] = log(gv.gvar('0.90(32)')) # excited omega
    prior['log(dE:rho)'][3] = log(gv.gvar('0.023(30)')) #excited rho

    prior.add('a:rhoo', [gv.gvar(0,0.1),gv.gvar(0,0.1),gv.gvar(0,0.1),gv.gvar(0,0.1)])
    prior.add('b:rhoo', [-prior['a:rhoo'][0],prior['a:rhoo'][1],prior['a:rhoo'][2],-prior['a:rhoo'][3]])
    prior.add('log(dE:rhoo)',[log(gv.gvar(1,1)) for i in range(nexp)])
    prior['log(dE:rhoo)'][0] = log(gv.gvar(2,2))
    ##

    return prior
##

def build_models(tmin,length,tag):
    """ build models """
    tdata = range(length)
    tfit = range(tmin,length-tmin) # all ts

    models = [

        Corr2(datatag=tag,tdata=tdata,tfit=tfit,tp=length,
            a=('a:rho','a:rhoo'),b=('b:rho','b:rhoo'),
            dE=('dE:rho','dE:rhoo'),s=(1,-1)),

    ]

    return models




periodic = False

print()
print('***** very coarse *****')

w0 = gv.gvar('0.1715(9)')
a_vcoarse = (w0/gv.gvar('1.13215(35)'))/0.197326968

# load very coarse correlator
dset = (gv.dataset.Dataset(data_dir + '/rho_vcphys_bothcharges_m0.001524.gpl'))

dsetavg = gv.dataset.avg_data(dset)


def single_fit(tag_, tstar) :

  dataf = dsetavg[tag_]
  tmin = 2

   # replace data after t* with fit
   #  t* can be changed


  fitter = CorrFitter(models=build_models(tmin,48,tag_))
  svdcut = 1e-6

  for nexp in [4]:
        print('=== nexp =',nexp,' tmin = ',tmin,' svdcut = ',svdcut)
        fit = fitter.lsqfit(data=dsetavg,prior=build_prior(nexp),p0=pfile,maxit=10000,svdcut=svdcut)#,fitter='gsl_multifit')
        print(fit)

  newdata = dataf[:tstar]
  #print("Before append/concatenate")
  #print(newdata, len(newdata))
  fitdata = Corr2(datatag=tag_,tdata=range(48),a=('a:rho','a:rhoo'),b=('b:rho','b:rhoo'),dE=('dE:rho','dE:rhoo'),s=(1,-1)).fitfcn(fit.p)
  newdata = np.append(np.array(newdata),fitdata[tstar:])
  #print(fitdata[tstar:], fitdata[tstar:].size)
  #print("After append/concatenate")
  #print(newdata, newdata.size)
  #sys.exit(0)

  # use Z_V from appendix of 1909.00756; divided by u_0 to make the 1-link currents used match
  vpol = g2.fourier_vacpol(newdata,Z=gv.gvar('0.93516(16)')/gv.gvar('0.820192(14)'),ainv=1/a_vcoarse,periodic=False)

  a_mu = g2.a_mu(vpol,Q=1./3.)

  print('a_mu:',a_mu, "for" , tag_, "t* = ", tstar)

  return a_mu 


for i in range(9,10):
  a_charge  = single_fit('charged-up', i) 
  a_neutral = single_fit('uncharged-up', i) 

  print("a_mu[QCD+QED] = " , a_charge)
  print("a_mu[QCD]     = " , a_neutral)

  rat = a_charge / a_neutral

  print("a_mu[QCD+QED]/a_mu[QCD] = " , rat)





