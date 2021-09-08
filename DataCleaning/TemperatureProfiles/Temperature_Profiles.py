'''
Python Routine to wrap together Temperature Map Fitting given a binning scheme
What should you already have?
    - binning scheme (WVT perhaps)
    - blank sky background
    - Fits and Image file used for spectral extraction
'''
#-----------------------------------------------IMPORTS--------------------------------------------------#
import os
import sys
from read_input import read_input_file
from Fits import PrimeFitting
from Plots import plot_data
from diffuse_specextract_blank import main_extract
#from deproject_mod import deproj_final
#-----------------------------------------------READ IN--------------------------------------------------#
inputs = read_input_file(sys.argv[1])
base = inputs['base_dir']+inputs['Name']
num_bins = len(inputs['reg_files'])
#----------------------------------------------------SPECTRA---------------------------------------------------#
if inputs['extract_spectrum'].lower() == 'true':
    for reg_file in inputs['reg_files']:
        # For each region extract the spectra in each ObsID
        main_extract(base,base+'/regions',inputs['ObsIDs'],reg_file)

#-----------------------------------------------FIT SPECTRA-------------------------------------------#
if inputs['fit_spec'].lower() == 'true':
    # Deprojection
    PrimeFitting( inputs['base_dir']+'/'+inputs['Name'],inputs['ObsIDs'],inputs['source_file'],int(num_bins),inputs['redshift'],inputs['n_H'],inputs['Temp_Guess'],inputs['Temp_data'],inputs['multi'])
