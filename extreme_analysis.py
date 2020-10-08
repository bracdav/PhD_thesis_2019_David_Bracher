# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:48:11 2020

@author: david
"""
import pandas as pd
import numpy as np
from glob import glob
from scipy.optimize import curve_fit
from thesis_db2019 import *
import thesis_db2019

path = 'D:\\Beamtime_SLS_raw\\bt1803 extreme\\CoO NPs\\'
folders = glob(path+'CoO*xmcd')

all_data = pd.DataFrame()

for folder in folders:
    files = glob(folder)  
    cols = ['#Ecrbk', 'FieldX', 'FieldZ', 'Pol', 'Temperature', 'NORMtey']                           
    data_extreme = XASExtreme(folder)
    data_extreme.list_files('*.txt')
    data_extreme.seperate_polarization_circ('\t')
    data_extreme.load_data(cols)
    xmcd = data_extreme.aggregate_data_circ()
    
    data_extreme.step_norm(interval = [0, 30, 470, 500])
    
    
    data_extreme.construct_intervall(ranges = [[0, 50], [51, 120], [121, 185], 
                                               [186, 250], [251, 499]],
                                     sample = 15)
    data_extreme.xas_match(int_given = True, intervall = data_extreme.intervall)
    data = data_extreme.construct_data()
    
    all_data = pd.concat((all_data, data), axis = 0)
    

all_data.sort_values(['temperature', 'fieldx'], ascending = [False, False], inplace = True)
all_data.loc[all_data.temperature < 30, 'temperature'] = 20
all_data.to_csv('extreme_data.csv')