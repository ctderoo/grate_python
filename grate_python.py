from numpy import *
import scipy.optimize as opt
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import csv
import pdb
import pickle
import os

'''
NOTE! The consoledir path must be correctly set to access the PCGrate solver, and should
be verified following installation.
'''
consoledir = 'C://Program Files (x86)//I. I. G. Inc//PCGrate-SX 6.1//ConsoleSolver//'
currentdir = os.getcwd()
    
def find_nearest_value(array,value):
    '''
    Utility function for closest lookup (used when interpolation isn't possible)
    Inputs:
    array - Searchable array
    value - value to be closely matched
    Outputs:
    array[ind] - the value in the array closest to the input 'value'.
    '''
    ind = abs(array-value).argmin()
    return array[ind]

def extract_calculation(filename,polarization = 'TE'):
    '''
    Function for reading and retrieving the efficiency information from the PCGrate
    calculation output file.So long as the inputXML file uses the NP polarization
    option, efficiencies for both polarizations will be calculated and a polarization
    can be specified. 
    '''
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
        
    return [wavelength,efficiency,order]

def convert_to_PCGrate_coords(graze,yaw):
    '''
    Converts the input coordinates often used in X-ray test geometry (graze,yaw) and converts
    them to the azimuth angle (phi) and polar angle (theta) used in the PCGrate coordinate
    system.
    Inputs:
    graze - incidence angle (in radians)
    yaw - yaw angle (in radians)
    Outputs:
    Phi - Azimuth angle for PCGrate (in radians)
    Theta - Polar angle for PCGrate (in radians)
    '''
    def equations(p,yaw,graze):
        [phi,theta] = p
        return (tan(yaw)-sin(theta)/tan(phi),sin(graze)-cos(theta)*cos(phi))
    
    solution = opt.fsolve(equations,[pi-graze,0.0],args=(yaw,graze))
    phi,theta = float(arcsin(sin(solution[0]))),float(arctan(tan(solution[1])))
    return phi,theta

def write_xml_file_lamellar(starting_file,phi,theta,period,depth,material_file = 'Au_CXRO_May_2006.ari',rms_roughness = 0.0,output_file = 'temp.xml'):
    '''
    Write a new xml file from a starting lamellar exemplar file. The starting exemplar file is set up initially
    in PCGrate.
    Inputs:
    starting_file - the exemplar xml file for the lamellar case. Note that this MUST be in the
    console solver directory.
    output_file - the desired path to the output xml file.
    phi - the PCGrate azimuthal coordinate to be written (output from convert_to_PCGrate_coords, should be in radians)
    theta - the PCGrate polar coordinate to be written (output from convert_to_PCGrate_coords, should be in radians)
    period - grating groove period (in nanometers)
    depth - total peak-to-valley height of the lamellar grating (in nanometers)
    material_file - input reflectivity file for the grating material
    Outputs:
    output_file - path to the written output xml file.
    ''' 
    tree = ET.parse(starting_file)
    root = tree.getroot()
    
    new_azimuth = str(phi*180/pi)
    new_polar = str(theta*180/pi)

    root[2][1].text = str(new_azimuth)
    root[2][2].text = str(new_polar)
    root[1][0].text = str(period)
    root[3][1][0][0][0].text = str(depth)
    root[3][2][2].text = str(rms_roughness)
    root[3][2][0].text = material_file
    temp_dir = os.getcwd()
    os.chdir(consoledir)
    tree.write(output_file,encoding = 'WINDOWS-1251')
    #tree.write(output_file)
    os.system(temp_dir)
    return output_file

def write_xml_file_blazed(starting_file,phi,theta,facet_angle,period,rms_roughness = 0.0,plateau_width = 5.0,output_file = 'temp.xml'):
    '''
    Write a new xml file from a starting blazed exemplar file. This currently does not have a
    working example, but is perhaps illustrative of how this works.
    '''
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
    temp_dir = os.getcwd()
    os.chdir(consoledir)
    tree.write(output_file)
    os.system(temp_dir)
    return output_file

def run_PCGrate_calculation(param_file,consoledir = consoledir):
    '''
    Function that actually runs the PCGrate calculation for the set-up input
    file. Returns the path to the results file for later reading by extract_calculation.
    Inputs:
    param_file - the parameter file that you'd like to have the calculation run on.
    Outputs:
    path - the path to the output calculation
    '''
    temp_dir = os.getcwd()
    os.chdir(consoledir)
    os.system('ConsoleSolver.exe -run')
    os.system('ConsoleSolver.exe -s ' + param_file + ' results.xml')
    while os.path.isfile(consoledir + 'results.xml') == False:
        pass
    os.chdir(temp_dir)
    return consoledir + '//results.xml'

def extract_order_diffraction_efficiency(fwave,feff,forder,order):
    '''
    Searches the output of extract_calculation for the efficiency/wavelength arrays of a
    specific order.
    Inputs:
    fwave - wavelength list from extract_calculation
    feff - efficiency list from extract_calculation
    forder - order list from extract_calculation
    order - the order you're searching for
    Outputs:
    wave,eff - the wavelength list, efficiency list of the matching order.
    '''
    ind = where(forder==order)
    eff = feff[ind]
    wave = fwave[ind]
    return wave,eff

def calculate_Littrow_eff(graze,facet_angle,groove_density,starting_file = 'blazed_grating_with_plateau.xml'):
    '''
    Calculates the efficiency of a grating put in the Littrow configuration.
    Inputs:
    graze - incidence angle (in degrees)
    facet_angle - the facet angle of the blazed grating. This facet angle is used to calculate
                the Littrow mounting for the grating.
    groove_density - the groove density of the blazed grating
    starting_file - the .xml file which will be used as a starting template to make the input .xml
                file for the calculation.
    '''
    groove_period = 1./groove_density*10**6
    grating_yaw = arcsin(tan(graze*pi/180)*tan(facet_angle*pi/180))
    
    os.chdir(consoledir)
    phi,theta = get_PCGrate_coordinates(graze*pi/180,grating_yaw)
    temp_file = write_xml_file(starting_file,'temp.xml',phi,theta,facet_angle,groove_period)
    results_path = run_PCGrate_calculation(temp_file)
    os.chdir(currentdir)
    return results_path

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