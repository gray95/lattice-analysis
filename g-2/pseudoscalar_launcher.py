
import numpy as np
import matplotlib.pyplot as plt
import math as m
from pseudoscalar_fit import fit_data 
import sys
sys.path.append('../code/plotters')
import gvar as gv
from plotlib import jackknife, calc_meff, plot_corr

from pseudoscalar_fit_params import OSC, CORRFIT, PLOT, SAVEFIG
from pseudoscalar_fit_params import corrpath, key, otherkey


if CORRFIT :
    fit_data(corrpath, key, otherkey)        # just a single channel fit

if PLOT :
    plot_corr(corrpath, SAVEFIG)

    

