import numpy as np
import os
import sys

filename_in = ['proc_corrs\hybrid\PION_PS_B_I_c_d_m0.8447\PION_PS_B_I_c_d_m0.8447.txt', 'proc_corrs\PION_PS_B_I_c_d_m0.8447b.txt', 'proc_corrs\PION_PS_B_I_c_d_m0.8447c.txt']
print(filename_in[0])
d = []
c = []

for i in range(3):
    print ("Reading data from  " , filename_in[i])
    d.append(open(filename_in[i], 'r'))
    c.append(np.genfromtxt(d[i], autostrip=True))

c[2] = c[2][5:]

#print(c[2])
    
end = np.add(c[0], (np.add(c[1],c[2])))*1/3

print(end.shape, end)

np.savetxt('Avg_PION_PS_B_I_c_d_m0.8447.txt', end)




