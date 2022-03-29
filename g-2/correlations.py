import gvar as gv
import sys
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

vc = '/home/gray/Desktop/lattice-analysis/data/qqed/vcoarse/Nml_rho_vc.gpl'
c3ml  = '/home/gray/Desktop/lattice-analysis/data/qqed/coarse/3ml_rho_coarse.gpl'

fnames = vc 
D  = gv.dataset.Dataset(fnames)
D_avg = gv.dataset.avg_data(D)


filepath = vc
f = open(filepath, 'r')
f1 = f.readlines()
ncfg = int(len(f1)/3)
ndel = 13

nbs = 50 # 
bs_test = (d for d in gv.dataset.bootstrap_iter(D,nbs))
for k in D.keys():
    print(next(bs_test)[k].shape)
bs_datalist = (gv.dataset.avg_data(d) for d in gv.dataset.bootstrap_iter(D,nbs))
corr = []
for B in bs_datalist:
    Y = {}
    for key in B:
        Y[key] = B[key] 
    cor = gv.evalcorr(Y)
    C = cor['VTvc3mlq0','VTvc3mlq1']
    corr.append(C)

print(len(corr))
corr_dist = gv.dataset.avg_data(corr)
org_corr_dist = gv.evalcorr(D_avg)

A = org_corr_dist['VTvc3mlq0','VTvc3mlq1'].diagonal()
B = corr_dist.diagonal()

print(B/A)

sys.exit(0)
#bootstrap_ensemble = []
#bootstrap_sample = gv.dataset.Dataset()

#Nsamples = range(ncfg)


#ax = sns.heatmap(C)


#plt.show()
#print(corr['VTvc3mlq0','VTc3mlq1'].diagonal())
#autocorr = gv.dataset.autocorr(D)

#avg = gv.dataset.avg_data(D, warn=True)
#corr = gv.evalcorr(avg)

