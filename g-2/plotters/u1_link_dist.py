

## This script plots a histogram of the average U(1) link values for selected cfgs

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# 50 cfgs (each) sampled randomly from the g-2 runs we have done
dan_u1_links = np.load('../misc/dan_u1_link_coarse.npy')
gray_u1_links = np.load('../misc/gray_u1_link_coarse.npy')

dlinks = [x.real for x in dan_u1_links]
glinks = [x.real for x in gray_u1_links]

dlinkmean = np.mean(dlinks)
glinkmean = np.mean(glinks)

dlinksdev = np.std(dlinks)
glinksdev = np.std(glinks)

diff = abs( (dlinkmean-glinkmean)/glinksdev )
s="diff = "+str(round(diff,4))+"$\sigma$"

plt.hist( dlinks, bins=50 , alpha=0.5, label="dans U(1) fields", color='red')
plt.hist( glinks, bins=50 , alpha=0.5, label="gauravs U(1) fields", color='blue')
plt.axvline(x=dlinkmean, color='red', linestyle='--', label='dan mean')
plt.axvline(x=glinkmean, color='blue', linestyle='--', label='gaurav mean')
plt.text(glinkmean+glinksdev, 4, s, fontsize=12)
plt.title('avg U(1) link on 50 randomly sampled cfgs')
plt.xlabel('link value')
plt.legend()


plt.tight_layout()
plt.savefig('../figures/u1_hist_link.png', dpi=500)
plt.show()
