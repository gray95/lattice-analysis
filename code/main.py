
import numpy as np
import matplotlib.pyplot as plt
import math as m
from corrfit import fit_data 
from plotlib import jackknife, calc_meff, plot_corr
import sys
import gvar as gv

from parameters import OSC, CORRFIT, PLOT, SAVEFIG
from parameters import corrpath, key, otherkey


if CORRFIT :
    fit_data(corrpath, key, otherkey)        # just a single channel fit

if PLOT :
    plot_corr(corrpath, key, SAVEFIG)

    

