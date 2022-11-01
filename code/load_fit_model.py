from __future__ import print_function   # makes this work for python2 and 3

import os
import sys
import dill as pickle
import bz2
#---------------------#
import collections
import gvar as gv
import numpy as np
#import corrfitter as cf
import datetime
from parameters import corrpath, e_str
from parameters import OSC, NOISE, WRITE_LOG
from parameters import SRCs, KEYFMT, tag, otag, SAVE_RES, save_res_name, save_fit
import matplotlib.pyplot as plt

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
plt.rcParams['savefig.dpi'] = 400

def main():
    print("Correlators - %s\nFit model - %s"%(corrpath,save_fit))
    fit = pickle.load(bz2.BZ2File(save_fit, 'rb'))

    print("reduced chi-squared: ", fit.chi2/fit.dof)

    fit.show_plots(save='../figures/'+e_str+'/{}.png', view='diff')

    sys.exit(0)
    # Error analysis 
    load_res = gv.load(save_res_name)
    fitp = load_res['fit.p']
    prior = load_res['prior']
    data = load_res['data']
    corrections = load_res['svdcut']

    E = np.cumsum(fitp[tag+'dE'])
    Eo = np.cumsum(fitp[otag+'dE'])
    print(gv.corr(E[0],Eo[0]))

    outputs = gv.load('./out/FINAL_RES')
#    outputs = collections.OrderedDict()    
#    outputs[tag+'grd'] = E[0]
#    outputs[otag+'grd'] = Eo[0]
    inputs=collections.OrderedDict()
    inputs['fitp'] = fitp
    inputs['prior'] = prior
    inputs['data']  = data
    inputs['svdcut'] = corrections
    print(gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs, colwidth=18, verify=True))


################################ FUNCTIONS ###########################################

        
if __name__ == '__main__':
    gv.ranseed(1234)
    if True:
        main()
    else:
        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()
        main()
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print (s.getvalue())
