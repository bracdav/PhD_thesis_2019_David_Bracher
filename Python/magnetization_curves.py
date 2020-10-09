# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:22:34 2020

@author: david
"""
from thesis_db2019 import *

filename = 'D:\\Beamtime_SLS\\python_scripts_final\\magnetic_moments.csv'

magn_curves = MagnetizationCurves(filename)
magn_curves.load_data()
magn_curves.set_zero()
magn_curves.plot_total(tp = 'total')
magn_curves.plot_total(tp = 'spin')
magn_curves.plot_total(tp = 'orbital')
magn_curves.fit_langevin()
magn_curves.construct_fit_lines()

plt.figure()
magn_curves.plot_total(tp = 'total')
magn_curves.plot_sims()
plt.xlabel('magnetic field (T)')
plt.ylabel('magnetic moments per atom (mu_b)')
plt.show()

np = MomentsPerNP(85)
np.set_volume()
np.set_unit_cell(4.2615)
np.set_interface(thickness = 0.42615*2)
np.get_moments(magn_curves.mag300[0])