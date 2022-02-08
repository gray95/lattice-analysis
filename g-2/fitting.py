

## the functions needed to fit the correlators

import os
import sys
import lsqfit
from corrfitter import Corr2,Corr3,CorrFitter
import numpy as np
from gvar import log,exp,evalcov
import gvar as gv
import g2tools as g2
import lsqfit as lsq

def build_prior(nexp, dp, a):
    prior = gv.BufferDict()

    prior.add('log(a1:vec:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(a1:vec:u)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(ao:vec:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(ao:vec:u)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(dE:vec:u)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:u)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:u)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:qed:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(a1:vec:qed:u)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(ao:vec:qed:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(ao:vec:qed:u)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(dE:vec:qed:u)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:qed:u)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:qed:u)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:qed:u)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(a1:vec:d)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(ao:vec:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(ao:vec:d)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(dE:vec:d)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:d)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:d)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:qed:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(a1:vec:qed:d)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(ao:vec:qed:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(ao:vec:qed:d)'][0] = log(gv.gvar(0.5*a,0.5*a*dp))
    prior.add('log(dE:vec:qed:d)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:qed:d)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:qed:d)',[log(gv.gvar(a,a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:qed:d)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(a1:vec:s)'][0] = log(gv.gvar(0.5, 0.5*dp))
    prior.add('log(ao:vec:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(ao:vec:s)'][0] = log(gv.gvar(0.5, 0.5*dp))
    prior.add('log(dE:vec:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:s)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:s)'][0] = log(gv.gvar(a,a*dp))

    prior.add('log(a1:vec:qed:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(a1:vec:qed:s)'][0] = log(gv.gvar(0.5,0.5*dp))
    prior.add('log(ao:vec:qed:s)',[log(gv.gvar(1,dp)) for i in range(nexp)])
    prior['log(ao:vec:qed:s)'][0] = log(gv.gvar(0.5,0.5*dp))
    prior.add('log(dE:vec:qed:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dE:vec:qed:s)'][0] = log(gv.gvar(a,a*dp))
    prior.add('log(dEo:vec:qed:s)',[log(gv.gvar(2*a,2*a*dp)) for i in range(nexp)])
    prior['log(dEo:vec:qed:s)'][0] = log(gv.gvar(a,a*dp))

    return prior

def fitargs(z):
    dp = z['p0']
    nexp = z['nexp']		
    prior = build_prior(nexp, dp)    
    return dict(prior=prior, fcn=model.fitfcn, data=ccdata[2:49])

def build_models(tag01,tag02,tmin,T):
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    tp = T # periodic
    models = [
        Corr2(datatag=tag01,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:s','ao:vec:s'),b=('a1:vec:s','ao:vec:s'),
            dE=('dE:vec:s','dEo:vec:s'),s=(1.,-1.)), 
        Corr2(datatag=tag02,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:s','ao:vec:qed:s'),b=('a1:vec:qed:s','ao:vec:qed:s'),
            dE=('dE:vec:qed:s','dEo:vec:qed:s'),s=(1.,-1.)),
    ]
    return models[0], models[1]

def build_models_ml(tag01,tag02,tag03,tmin,T):
 
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    
    tp = T # periodic
    models = [
        Corr2(datatag=tag01,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),
            dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)),
        Corr2(datatag=tag02,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),
            dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)),
        Corr2(datatag=tag03,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),
            dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.))
    ]
    return models

def build_models_phys(tag01,tag02,tag03,tag04,tmin,T):
 
    tdata = range(T)
    tfit = range(tmin,T+1-tmin) # all ts
    
    tp = T # periodic
    models = [
        Corr2(datatag=tag01,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:u','ao:vec:u'),b=('a1:vec:u','ao:vec:u'),
            dE=('dE:vec:u','dEo:vec:u'),s=(1.,-1.)),
        Corr2(datatag=tag02,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:u','ao:vec:qed:u'),b=('a1:vec:qed:u','ao:vec:qed:u'),
            dE=('dE:vec:qed:u','dEo:vec:qed:u'),s=(1.,-1.)),
        Corr2(datatag=tag03,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:d','ao:vec:d'),b=('a1:vec:d','ao:vec:d'),
            dE=('dE:vec:d','dEo:vec:d'),s=(1.,-1.)),
        Corr2(datatag=tag04,tp=tp,tdata=tdata,tfit=tfit,
            a=('a1:vec:qed:d','ao:vec:qed:d'),b=('a1:vec:qed:d','ao:vec:qed:d'),
            dE=('dE:vec:qed:d','dEo:vec:qed:d'),s=(1.,-1.))
    ]
    return models



def make_data(filename,norm=1):
    dset = gv.dataset.Dataset(filename)
    dset = gv.dataset.bin_data(dset, binsize=1)
    s = gv.dataset.svd_diagnosis(dset)
    print(dset.keys())
    for tag in dset.keys():
        dset[tag] = norm*np.array(dset[tag])
        #dset[tag] = dset[tag][:3]
        print("number of cfgs| ", tag, " : ", dset[tag].shape[0]) 
    cdata = gv.dataset.avg_data(dset)
    return (cdata, dset.keys(), s.svdcut)

def print_results(fit, ainv, Zv, Zvqed):
    p = fit.p
    E = np.cumsum(p['dE:vec:s'])
    a = p['a1:vec:s']
    Eo = np.cumsum(p['dEo:vec:s'])
    ao = p['ao:vec:s']

    Eqed =  np.cumsum(p['dE:vec:qed:s'])
    aqed = p['a1:vec:qed:s']

    inputs = {'E[0]':E[0], 'Eqed[0]':Eqed[0], 'a[0]':a[0], 'aqed[0]':aqed[0], 'a':1/ainv, 'Zv':Zv, 'Zvqed':Zvqed }
    phi_mass = E[0] * ainv
    phi_qed_mass = Eqed[0] * ainv
    phi_mass_diff = (phi_qed_mass-phi_mass)*1000
    phi_mass_rt = phi_qed_mass/phi_mass
    phi_mass_diff_rt = phi_mass_diff/phi_mass
    print("phi mass is %s GeV in qcd"%(phi_mass))
    print("phi mass is %s GeV in qcd+qqed"%(phi_qed_mass))
    print("phi mass diff = %s MeV"%((phi_mass_diff)))
    print("phi mass ratio = %s\n"%(phi_mass_rt))

    print("a[0] is %s\n"%(a[0]))

    f_phi_lat = a[0] * gv.sqrt( 2.0/E[0] )
    f_phi = 1000 * f_phi_lat * ainv * Zv 
    f_phi_qed_lat = aqed[0] * gv.sqrt( 2.0/Eqed[0] )
    f_phi_qed = 1000 * f_phi_qed_lat * ainv * Zvqed
    f_rt = f_phi_qed/f_phi 
    f_diff = f_phi_qed - f_phi 

    print("vector leptonic decay constant is %s in lattice units "%(f_phi_lat)) 
    print("vector leptonic decay constant is %s in lattice units "%(f_phi_qed_lat)) 
    print("vector leptonic decay constant is %s in MeV "%(f_phi)) 
    print("vector leptonic decay constant is %s in MeV "%(f_phi_qed)) 
    print("vector leptonic decay diff is %s MeV"%(f_diff)) 
    print("vector leptonic decay constant ratio is %s\n"%(f_rt)) 

    outputs = { 'phi M':phi_qed_mass, 'phi f':f_phi_qed, 'M/M0':phi_mass_rt, 'f/f0':f_rt, 'M-M0':phi_mass_diff, 'f-f0':f_diff}
    print(gv.fmt_errorbudget(outputs=outputs, inputs=inputs))
    
