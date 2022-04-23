import gvar as gv
import sys
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

vc = '/home/gray/Desktop/lattice-analysis/data/qqed/vcoarse/vt_vcoarse.gpl'
c  = '/home/gray/Desktop/lattice-analysis/data/qqed/coarse/vt_coarse.gpl'

TAG1 = 'VTc7mlq2'
TAG2 = 'VTc3mlq0'
TAGS = ['VTvc3mlq0', 'VTvc5mlq0', 'VTvc7mlq0', 'VTvcmsq0']

fnames = c 
D  = gv.dataset.Dataset(fnames)
D_avg = gv.dataset.avg_data(D)
org_corr_dist = gv.evalcorr(D_avg)
nkeys = len(list(D.keys()))
ran_key = list(D.keys())[3]
T = D_avg[ran_key].shape[0]

f = open(fnames, 'r')
f1 = f.readlines()
ncfg = int(len(f1)/nkeys)
print("%d cfgs in %s"%(ncfg,fnames))

nbs = 50
print('%d bootstrap samples'%nbs)
bs_datalist = (gv.dataset.avg_data(d) for d in gv.dataset.bootstrap_iter(D,nbs))
bs_corr = gv.dataset.Dataset()
for B in bs_datalist:
    cor = gv.evalcorr(B)
    for K in cor:
        bs_corr.append(K, cor[K].diagonal())

bs_corr_avg = gv.dataset.avg_data(bs_corr, bstrap=True)

A = org_corr_dist[TAG1,TAG2].diagonal()
B = bs_corr_avg[TAG1, TAG2]


######################## PLOTTING ###########################################
plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('axes', linewidth=0.5)


#plt.title('correlation of '+TAG1+' with '+TAG2)
#plt.errorbar(range(T), [y.mean for y in B], xerr=0, yerr=[y.sdev for y in B], fmt='h', color='blue', alpha=0.5, label='bootstrap correlation, nbs='+str(nbs))
#plt.errorbar(range(T), A, xerr=0, yerr=0, fmt='h', color='red', alpha=0.5, label='total sample')

#handles,labels = plt.gca().get_legend_handles_labels()
#handles = [h[0] for h in handles]
#plt.legend(handles=handles,labels=labels,frameon=False,loc='lower right')


sns.heatmap(org_corr_dist[TAG1,TAG2][::2,::2])
plt.show()
#autocorr = gv.dataset.autocorr(D)
plt.tight_layout()
#avg = gv.dataset.avg_data(D, warn=True)
#corr = gv.evalcorr(avg)

