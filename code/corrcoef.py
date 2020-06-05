

# compute how closely two different correlators are correlated

import numpy as np
import os 
import sys

file1 = os.path.join('proc_corrs','l3264f211b600m00507m0507m628a-coul-v5','Hy_m0.650_p000.txt')
file2 = os.path.join('proc_corrs','l3264f211b600m00507m0507m628a-coul-v5','t32_Hy_m0.650_p000.txt')

corr1 = np.genfromtxt(file1, max_rows=5, usecols=range(1,65))
corr2 = np.genfromtxt(file2, max_rows=5, usecols=range(1,65))

#corr1 = corr1[0,1]
#corr2 = corr2[0,1]

cov = np.corrcoef(corr1, corr2, rowvar=False)
covariance_matrix = np.cov(corr1, corr2, rowvar=False)

print(cov)



