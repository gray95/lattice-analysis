To extract the data from the correlates in the 
subdirectory of work you need to modify the 
python script: reformat.py

The corr_tag variable needs to be updated to get the
required channel:

corr_tag =  'PION_PS_C_C_d_d_m0.838_m0.838_p000'
corr_tag =  'PION_SC_C_C_d_d_m0.838_m0.838_p000'
corr_tag =  'RHO_VT_C_C_d_d_m0.838_m0.838_p000'
corr_tag =  'RHO_PV_C_C_d_d_m0.838_m0.838_p000'


The command 
  print ("1-+_4ape " , end = "")
also should be updated. 
