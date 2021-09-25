
import numpy as np
import matplotlib.pyplot as plt
import math as m
from corrfit import fit_data 
import sys
sys.path.append('plotters')
import gvar as gv
from plotlib import jackknife, calc_meff, plot_corr

from parameters import OSC, CORRFIT, PLOT, SAVEFIG
from parameters import corrpath, key, otherkey


if CORRFIT :
    fit_data(corrpath, key, otherkey)        # just a single channel fit

if PLOT :
    plot_corr(corrpath, SAVEFIG)

    

