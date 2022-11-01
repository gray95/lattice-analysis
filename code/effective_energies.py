import gvar as gv
import numpy as np
from commonweal import onemp_paths
import corrfitter as cf
import matplotlib.pyplot as plt
import sys

def smoothed_energies( E ):
    smooth = [ 0.25*(E[t+1]-2*E[t]+E[t-1]) for t in range(1,len(E)-1) ]
    return smooth


e_str = 'vc'
SIZE  = {'vc':48, 'c':64, 'f':96}

T = SIZE[e_str]

corr = gv.dataset.avg_data(cf.read_dataset(onemp_paths[e_str]))

C = np.array( [[[corr['onemp.ll'][t],corr['onemp.lg'][t]], [corr['onemp.gl'][t],corr['onemp.gg'][t]]] for t in range(T)] )

eigval = []

print("eigenvalues")
for t in range(T):
    eigval.append(gv.linalg.eigh(C[t])[0][-1])

sort_eigval = np.sort(eigval)

NO_E = [gv.log((eigval[t+1]/eigval[t])) for t in range(20)] # Non Oscillating effective energies 
print(NO_E)

print(smoothed_energies(NO_E))

sys.exit(0)

##################################################################################
################## PLOT ################################

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=0.5)
#plt.rcParams['savefig.dpi'] = 200
plt.figure(figsize=((4.5,4.5/1.618)))
plt.gca().tick_params(right=True,top=True,direction='in')


plt.errorbar(x=range(T), y=[y.mean for y in sort_eigval],yerr=[y.sdev for y in sort_eigval], fmt='h')
#plt.yscale('log')
plt.show()
