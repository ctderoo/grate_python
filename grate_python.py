from numpy import *
import scipy.optimize as opt
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import csv
import pdb
import pickle
import os

consoledir = 'C://Program Files (x86)//I. I. G. Inc//PCGrate-SX 6.1//ConsoleSolver//'
currentdir = os.getcwd()
    
def find_nearest_value(array,value):
    ind = abs(array-value).argmin()
    return array[ind]

def scanwave(filename,polarization = 'TE'):
    tree = ET.parse(filename)
    root = tree.getroot()

    wavelength = []
    TM_efficiency = []
    TE_efficiency = []
    order = []
    for step in root.iter('Step'):
        for ordloop in step.iter('Order'):
            wavelength.append(float(step.get('Scanning_parameter_value')))
            order.append(int(ordloop.get('Order_number')))
            TM_efficiency.append(float(ordloop.get('Efficiency_TM')))
            TE_efficiency.append(float(ordloop.get('Efficiency_TE')))
            
    wavelength = array(wavelength)
    order = array(order)
    TM_efficiency = array(TM_efficiency)
    TE_efficiency = array(TE_efficiency)
    
    if polarization == 'TE':
        efficiency = TE_efficiency
    if polarization == 'TM':
        efficiency = TM_efficiency
    if polarization == 'NP':
        efficiency = (TE_efficiency+TM_efficiency)/2
        
    return [wavelength,order,efficiency]

def get_PCGrate_coordinates(gamma,grating_yaw):
    def equations(p,yaw,gamma):
        [phi,theta] = p
        return (tan(yaw)-sin(theta)/tan(phi),sin(gamma)-cos(theta)*cos(phi))
    
    solution = opt.fsolve(equations,[pi-gamma,0.2],args=(grating_yaw,gamma))
    phi,theta = float(arcsin(sin(solution[0]))),float(arctan(tan(solution[1])))
    return phi,theta

def write_xml_file(starting_file,output_file,phi,theta,facet_angle,period,rms_roughness = 0.0,plateau_width = 5.0):
    tree = ET.parse(starting_file)
    root = tree.getroot()
    
    new_azimuth = str(phi*180/pi)
    new_polar = str(theta*180/pi)

    root[2][1].text = str(new_azimuth)
    root[2][2].text = str(new_polar)
    root[1][0].text = str(period)
    root[3][1][0][0][0].text = str((period - plateau_width)/(1./tan(facet_angle*pi/180) + 1./tan((180 - 70.5 - facet_angle)*pi/180)))            # If you ever want to scan in the plateau width space.
    root[3][1][0][0][1].text = str(facet_angle)                 # If you ever want to scan in facet angle space.
    root[3][1][0][0][2].text = str(180 - 70.6 - facet_angle)    # To get the blazed silicon grating profile correct.
    root[3][2][2].text = str(rms_roughness)
    tree.write(output_file)
    return output_file

def run_PCGrate_calculation(param_file,consoledir = consoledir):
    os.system('ConsoleSolver.exe -run')
    os.system('ConsoleSolver.exe -s ' + param_file + ' ' + param_file + '_results.xml')
    while os.path.isfile(consoledir + param_file + '_results.xml') == False:
        pass
    return consoledir + param_file + '_results.xml'

def calculate_Littrow_eff(gamma,facet_angle,groove_density,starting_file = 'blazed_grating_with_plateau.xml'):
    '''
    Calculates the efficiency of a grating put in the Littrow configuration.
    Inputs:
    gamma - incidence angle (in degrees)
    facet_angle - the facet angle of the blazed grating. This facet angle is used to calculate
                the Littrow mounting for the grating.
    groove_density - the groove density of the blazed grating
    starting_file - the .xml file which will be used as a starting template to make the input .xml
                file for the calculation.
    '''
    groove_period = 1./groove_density*10**6
    grating_yaw = arcsin(tan(gamma*pi/180)*tan(facet_angle*pi/180))
    
    os.chdir(consoledir)
    phi,theta = get_PCGrate_coordinates(gamma*pi/180,grating_yaw)
    temp_file = write_xml_file(starting_file,'temp.xml',phi,theta,facet_angle,groove_period)
    results_path = run_PCGrate_calculation(temp_file)
    os.chdir(currentdir)
    return results_path

def extract_order_diffraction_efficiency(fwave,feff,forder,order):
    ind = where(forder==order)
    eff = feff[ind]
    wave = fwave[ind]
    return wave,eff

def plot_diffraction_efficiency(fwave,feff,forder,order_range=None):
    fig = plt.figure(figsize = (16,12))
    if order_range == None:
        order_range = arange(amin(forder),amax(forder)+1)
    for order in order_range:
        wave,eff = extract_order_diffraction_efficiency(fwave,feff,forder,order)
        plt.plot(1240./wave,eff,label = 'Order ' + str(order))
    plt.ylabel('Efficiency')
    plt.xlabel('Energy (eV)')
    plt.legend()
    return fig,plt.gca()

# def plot_sum_of_orders_efficiency(fwave,feff,forder,energy_range = linspace(300,1300,(1300-300)/50 + 1)):
#     sum_of_orders = []
#     for energy in energy_range:
#         wave = 1240./energy             # Assuming energy is specified in eV, wave in nm
#         wave_ind = argmin(abs(wave-fwave))
#         lookup_wave = fwave[0][wave_ind]
#         sum_of_orders.append(sum(feff[where(fwave == lookup_wave)]))
#     plt.plot(energy_range,asarray(sum_of_orders),'b--',marker = 'o',markerfacecolor = 'b',markersize = 10)
#     return sum_of_orders,energy_range

########################################################
# def identify_closest_diff_eff_calc(gamma,facet_angle,groove_density):
#     all_library_files = os.listdir('C:\Users\Casey DeRoo\Software\Diffraction_Efficiency_Libraries')
#     all_library_files = [s for s in all_library_files if 'gamma' in s]                                                                  # Scrub the list of non-library files
#     
#     calculated_gammas = list(set([filename[filename.index('gamma_')+6:filename.index('_facet_angle')] for filename in all_library_files]))
#     chosen_gamma = find_nearest_value(array(map(float,calculated_gammas)),gamma)
#     acceptable_files = [s for s in all_library_files if str(chosen_gamma) in s]
#     
#     calculated_facet_angles = list(set([filename[filename.index('facet_angle_')+12:filename.index('_grden')] for filename in acceptable_files]))
#     chosen_facet_angle = find_nearest_value(array(map(float,calculated_facet_angles)),facet_angle)
#     acceptable_files = [s for s in acceptable_files if str(chosen_facet_angle) in s]
#     
#     calculated_grden = list(set([filename[filename.index('_grden_')+7:filename.index('_Littrow')] for filename in acceptable_files]))
#     chosen_grden = find_nearest_value(array(map(float,calculated_grden)),groove_density)
#     acceptable_files = [s for s in acceptable_files if str(chosen_grden) in s]
#     
#     acceptable_file = acceptable_files[0]
#     return acceptable_files
# Identifying the correct theoretical file from the directory.
# material_name = 'SiO2'
# os.chdir(theo_directory)
# right_file = os.listdir('.')
# right_file = [f for f in right_file if 'materialSiO2' in f]
# #right_file = [f for f in right_file if 'grden_6000' in f]
# #right_file = [f for f in right_file if 'gamma_1.5' in f]
# #right_file = [f for f in right_file if 'facet_angle_54' in f]
# right_file = [f for f in right_file if 'TE_Littrow.npy' in f]
# 
# [fwave,feff,forder] = load(right_file[0])
# os.chdir(home_directory)
# 
# # Make theoretical efficiencies by orders plot.
# plt.ion()
# plt.figure(figsize=(16,10))
# plot_diffraction_efficiency(1240./fwave,feff,forder,order_range = arange(-22,4,1))
# plt.legend(fontsize = 10,ncol=2)
# plt.title('Predicted Diffraction Efficiencies for\nPitch = 160 nm, ' + r'$g$' + ' = 1.5 deg, ' + r'$\delta$' + ' = 54.7 deg', fontsize = 24)
# plt.xlabel('Energy (eV)',fontsize = 20)
# plt.ylabel('Efficiency', fontsize = 20)
# plt.savefig(material_name + '_Theoretical_Efficiencies_ByOrder_ExpandedEnergyRange.png')
# plt.close()
# 
# plt.ion()
# plt.figure(figsize=(16,10))
# plot_diffraction_efficiency(1240./fwave,feff,forder,order_range = arange(-22,4,1))
# 
# handles, labels = plt.gca().get_legend_handles_labels()
# newLabels, newHandles = ['-13th','-12th','-11th','-10th','-9th','-8th','-7th','-6th','-5th','-4th','-3rd','-2nd','-1st','0th','+1st'], []
# 
# for handle, label in zip(handles, labels):
#     if label in newLabels:
#         newHandles.append(handle)
#     #
#     #if label not in newLabels:
#     #    newLabels.append(label)
#     #    newHandles.append(handle)
# plt.legend(newHandles, newLabels,fontsize=10,ncol = 2)
# 
# #plt.legend(fontsize = 10)
# plt.title('Predicted Diffraction Efficiencies for\nPitch = 160 nm, ' + r'$g$' + ' = 1.5 deg, ' + r'$\delta$' + ' = 54.7 deg', fontsize = 24)
# plt.xlabel('Energy (eV)',fontsize = 20)
# plt.ylabel('Efficiency', fontsize = 20)
# plt.xlim(300,1300)
# plt.savefig(material_name + '_Theoretical_Efficiencies_ByOrder_ALSRange.png')
# plt.close()
# 
# # Make 'Sum-of-Orders' efficiency plot.
# plt.ion()
# plt.figure(figsize = (16,10))
# sum_of_orders,energy_range = plot_sum_of_orders_efficiency(fwave,feff,forder)
# plt.title('Predicted Summed Efficiency for\nPitch = 160 nm, ' + r'$g$' + ' = 1.5 deg, ' + r'$\delta$' + ' = 54.7 deg', fontsize = 24)
# plt.xlabel('Energy (eV)',fontsize = 20)
# plt.ylabel('Summed Efficiency', fontsize = 20)
# plt.xlim(300,1300)
# plt.savefig(material_name + '_Theoretical_Efficiencies_Sum.png')
# plt.close()
# 
# 
# 
# 
