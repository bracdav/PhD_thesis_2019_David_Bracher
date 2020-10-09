# PhD_thesis_2019_David_Bracher

#Disclaimer:
#This codebase was used in order to analyse data used for my PhD thesis (
#'Antiferromagnetic properties of 3d metaloxide nanoparticles' - David Bracher,
#University of Basel/Paul Scherrer Institut, 2019). The PhD thesis will be
#available public by  the library of the University of Basel:
#https://edoc.unibas.ch/view/type/doctoral.date.html

#The codebase consists of a library containing all essential functions and objects. It contains code to analyse x-ray absorption spectra with linear polarisation and x-ray magnetic linear dichroism spectra recorded by x-ray photoemission electron microscopy (at the Surface/Interface Microscopy (SIM) beamline at the Swiss Light Source at Paul Scherrer Institut) as well as non-spatially resolved XAS and x-ray magnetic circular dichroism spectra (obtained at Extreme beamline at the Swiss Light Source).

#At the moment the codes to analyse the x-ray diffraction and reflection high energy electron diffraction data are not contained in the code base. However, those codes will also be translated from their initial Matlab versions soon. Although the Matlab codes are included for now, which will be obsolete when  the Python codes are finished.

#The library is called thesis_db2019 contains:
#The library should be imported calling 'from thesis_db2019 import *' in order to be able to call the calsses and functions directly. Further, importing the library this way will directly import the folowing packages:
	#glob from glob
	#matplotlib.pyplot as plt
	#numpy as np
	#pandas as pd
	#curve_fit from scipy.optimize
import scipy.constants as nat
import seaborn as sns
from Stoner.Fit import langevin as lv

#Classes:
	
	#CoO:
		#Inherits all functionality from Class DichroicSpectrum. Specific 		functionality for analysis of CoO nanoparticles can be added here.

	#DichroicSpectrum:
		#calc_dichroism
		#calc_dichroism_norm
		#norm_on_feature
		#select_feature
		#xas_match
	
	#FeOOH:
		#construct_fit
		#det_features
		#execute_fit
		#fit_features_l2
		#fit_features_l23
		#fit_features_l3
		#load_files
		#sim_lines

	#MagnetizationCurves:
		#construct_fit_lines
		#fit_langevin
		#plot_sims
		#plot_total
		#set_zero

	#MomentsPerNP:
		#get_moments
		#set_interface
		#set_unit_cell
		#set_volume
		
	#SpinAxis:
		#construct_fit
		#execute_fit
		#load_files
		#sim_lines

	#SumRules:
		#calculate_A
		#calculate_B
		#calculate_C
		#calculate_orb
		#calculate_spin
		#get_param_pairs
		#use_parm

	#XASectra:
		#construct_nums
		#load_energy
		#xas_load
		#xas_mean
		#xas_normalize

	#XASExtreme:
		#aggregate_data_circ
		#construct_data
		#construct_interval
		#load_data
		#list_files
		#seperate_polarization
		#step_norm
		#xas_match