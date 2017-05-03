from numpy import *
import matplotlib.pyplot as plt
import os
import pdb

import grate_python.grate_python as gp

##################################################
# Grating parameters
##################################################
graze = 1.5*pi/180
yaw = 0.0*pi/180
period = 400.0 # in nanometers
depth = 50.0 # in nanometers
orders = arange(7) + 2  
##################################################
home_directory = os.getcwd()
exemplar_file = 'C:\\Program Files (x86)\\I. I. G. Inc\\PCGrate-SX 6.1\\ConsoleSolver\\lamellar_grating_conical.xml'
output_file = 'C:\\Program Files (x86)\\I. I. G. Inc\\PCGrate-SX 6.1\\ConsoleSolver\\temp.xml'

phi,theta = gp.convert_to_PCGrate_coords(graze,yaw)
param_filename = gp.write_xml_file_lamellar(exemplar_file,phi,theta,period,depth)
results_fn = gp.run_PCGrate_calculation(param_filename)
[wavelength,efficiency,order] = gp.extract_calculation(results_fn)
fig,ax1 = gp.plot_diffraction_efficiency(wavelength,efficiency,order,order_range = orders)

mgk_energy,mgk_eff = 1250,[]
ax1.vlines(mgk_energy,0,0.1,'k',linestyle = 'dashed')
for wanted_order in orders:
	fwave,feff = gp.extract_order_diffraction_efficiency(wavelength,efficiency,order,wanted_order)
	mgk_eff.append(feff[argmin(abs(fwave - 1240./mgk_energy))])
	



