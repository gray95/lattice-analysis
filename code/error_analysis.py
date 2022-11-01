import matplotlib
##matplotlib.use('Agg')
import numpy as np
import gvar as gv
import lsqfit as lsq
import matplotlib.pyplot as plt
import sys
sys.path.append('./plotters')
from results import onemp_m, pp_m
from commonweal import ainv, hbarc
from parameters import save_res_name
import collections


# error budget
#print("Correlators - %s\nFit model - %s"%(corrpath,save_fit))
#fit = pickle.load(bz2.BZ2File(save_fit, 'rb'))
#print("reduced chi-squared: ", fit.chi2/fit.dof)
#print("t0 is {}\nfit-range [{},{}]".format(t0, TFIT[0], TFIT[-1]))
load_res = gv.load(save_res_name)
print(load_res['fit.p'])
outputs = collections.OrderedDict()
outputs['cont onemp'] = continuum_limit_h[0]
outputs['cont pp'] = continuum_limit_pp[0]
inputs = collections.OrderedDict()
inputs['onemp data'] = onemp_am['vc']
inputs['prior'] = prior
print(gv.fmt_values(outputs))
print(gv.fmt_errorbudget(outputs, inputs, colwidth=18))


