# VARIOUS SMALL METHODS TO GET A FEEL FOR THE DATA
# compute how closely two different correlators are correlated
# compute the fractional error in the correlator

import gvar as gv
import numpy as np
import os 
import sys
import matplotlib.pyplot as plt

bp = '../data/proc_corrs'
conf1= 'l3296f211b630m0074m037m440-coul-v5'
conf2= 'l3296f211b630m0074m037m440-coul-v5'
#l3296f211b630m0074m037m440-coul-v5-new

filenames = ['PION_PS_m0.450.txt', 't0_onempHy_m0.450.txt']

file1 = os.path.join(bp, conf1, filenames[0])
file2 = os.path.join(bp, conf2, filenames[1])

corr1 = gv.dataset.Dataset(file1)
corr1 = corr1.toarray()
corr1 = corr1['pion.ll']

corr2 = gv.dataset.Dataset(file2)
corr2 = corr2.toarray()
corr2 = corr2['onemp.gg']

#cov = np.corrcoef(corr1, corr2, rowvar=False)
#covariance_matrix = np.cov(corr1, corr2, rowvar=False)
#print(cov)

corr1_unc = np.std(corr1, axis=0)
corr1_mean= np.mean(corr1, axis=0)

corr2_unc = np.std(corr2, axis=0)
corr2_mean= np.mean(corr2, axis=0)

frac_err = np.divide(corr1_unc, np.abs(corr1_mean))

print(frac_err)

plt.plot(frac_err)
plt.show()


