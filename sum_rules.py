# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:21:14 2020

@author: david
"""
from thesis_db2019 import *


filename = 'D:\\Beamtime_SLS\\python_scripts_final\\extreme_data.csv'

moments = SumRules(filename)
moments.load_data()
moments.get_param_pairs()

ind = [0,2,4,6]
cols = ['20K', '20K_err', '50K', '50K_err', '300K', '300K_err']
spin_df = pd.DataFrame(index = ind, columns = cols)
orb_df = pd.DataFrame(index = ind, columns = cols)
tot_df = pd.DataFrame(index = ind, columns = cols)
for temp, h in moments.params:
    moments.use_parm(temp = temp, h = h)
    moments.calculate_A([30, 150])
    moments.calculate_B([150, 300])
    moments.calculate_C([30, 300])
    spin = moments.calculate_spin()
    orb = moments.calculate_orb()
    tot = np.abs(spin + orb)
    moments.calculate_A([30, 150], off = 0.001)
    moments.calculate_B([150, 300], off = 0.001)
    moments.calculate_C([30, 300], off = 0.001)
    spin_err = spin - moments.calculate_spin()
    orb_err = orb - moments.calculate_orb()
    tot_err = np.abs(np.mean([spin_err, orb_err]))
    
    spin_df.loc[h,str(int(temp))+'K'] = np.abs(spin)
    spin_df.loc[h,str(int(temp))+'K_err'] = spin_err
    spin_df['type'] = 'spin_moment'
    
    orb_df.loc[h,str(int(temp))+'K'] = np.abs(orb)
    orb_df.loc[h,str(int(temp))+'K_err'] = orb_err
    orb_df['type'] = 'orbital_moment'
    
    tot_df.loc[h,str(int(temp))+'K'] = tot
    tot_df.loc[h,str(int(temp))+'K_err'] = tot_err
    tot_df['type'] = 'total_moment'
    
res_df = pd.concat((spin_df, orb_df, tot_df), axis = 0) 
res_df.reset_index(inplace = True)
cols = ['H', '20K', '20K_err', '50K', '50K_err', '300K', '300K_err', 'type']
res_df.columns = cols
res_df.set_index('H', inplace = True)
res_df.to_csv('magnetic_moments.csv')
print(res_df)