'''
------------------------------------------------------
GOAL:
    Step through bins (spectra) and calculate the temperature value
    of each bin using XSPEC
------------------------------------------------------
INPUTS:
    dir - Full Path to Main Directory (e.g. '/home/user/Documents/Chandra/12833/repro/binned/')
    file_name - FIlename of PI/PHA spectrum (e.g. 'imageA')
    output_file - Filename for output containing temperature information (e.g. 'Temp_bin')
    num_files - number of bins (e.g. 100)
    redshift - redshift of object (e.g. 0.659)
    n_H - Hydrogen Equivalent Column Density in units of 10^{22} atoms cm^{-2} (e.g. 3.6e-2)
    Temp_guess - Guess for Temperature value in units of KeV (e.g. 5)
------------------------------------------------------
OUTPUTS:
    A file which contains the bin number and associated
    temperature and reduced chi-squared
------------------------------------------------------
ADDITIONAL NOTES:
    Be sure to run heainit first
------------------------------------------------------
'''
#from astropy.io import fits
import os
import subprocess
from sherpa.astro.xspec import *
from sherpa.astro.all import *
from sherpa.astro.ui import *
from sherpa.all import *
cflux = XScflux()

#TURN OFF ON-SCREEN OUTPUT FROM SHERPA
import logging
logger = logging.getLogger("sherpa")
logger.setLevel(logging.WARN)
logger.setLevel(logging.ERROR)
def set_log_sherpa():
    p = get_data_plot_prefs()
    p["xlog"] = True
    p["ylog"] = True
    return None

def isFloat(string):
    if string == None:
        return False
    try:
        float(string)
        return True
    except ValueError:
        return False
    except TypeError:
        return False

#------------------------------INPUTS------------------------------------------#

#------------------------------------------------------------------------------#
#Dynamically set source for OBSID
def obsid_set(src_model_dict,bkg_model_dict,obsid,bkg_src, obs_count,redshift,nH_val,Temp_guess, spec_count, fit_params, AGN):
    '''
    Function to set the source and background model for the observation

    Args:
        src_model_dict : Dictionary of all the source models set for particular region
        bkg_model_dict : Dido except for the background models
        obsid : The source pha/pi file
        bkg_src : The background pha/pi file
        obs_count : The count of the current observation

    '''
    temperature, abundance = fit_params
    load_pha(obs_count,obsid,use_errors=True) #Read in
    if obs_count == 1:  # Set if this is the first
        if spec_count < AGN:
            src_model_dict[obsid] = xsphabs('abs' + str(obs_count)) * (xsapec(
                'apec' + str(obs_count)) + xszphabs('zphabs' + str(obs_count)) * (xspowerlaw('pow'+ str(obs_count)) + xsgaussian('gauss' + str(obs_count))))  # set model and name
            thaw(get_model_component('zphabs1').redshift)
            get_model_component('gauss1').LineE = 6.4
            get_model_component('pow1').PhoIndex = 1.8
        else:
            src_model_dict[obsid] = xsphabs('abs'+str(obs_count)) * xsapec('apec'+str(obs_count)) #set model and name
        # Change src model component values
        get_model_component('apec' + str(obs_count)).kT = temperature
        get_model_component('apec' + str(obs_count)).redshift = redshift  # need to tie all together
        get_model_component('apec' + str(obs_count)).Abundanc = abundance
        thaw(get_model_component('apec' + str(obs_count)).Abundanc)
        get_model_component('abs1').nH = nH_val  # change preset value
        freeze(get_model_component('abs1'))

    set_source(obs_count, src_model_dict[obsid]) #set model to source
    set_bkg(obs_count, unpack_pha(bkg_src))
    bkg_model_dict[obsid] = xsapec('bkgApec'+str(obs_count))+get_model_component('abs1')*xsbremss('brem'+str(obs_count))
    set_bkg_model(obs_count,bkg_model_dict[obsid])
    #Change bkg model component values
    get_model_component('bkgApec' + str(obs_count)).kT = 0.18
    freeze(get_model_component('bkgApec'+str(obs_count)).kT)
    get_model_component('brem' + str(obs_count)).kT = 40.0
    freeze(get_model_component('brem' + str(obs_count)).kT)

    return None
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def sherpa_model(fit_params, obs_ct, redshift, nH_val):
    temperature, abundance = fit_params
    mdl = xsphabs('abs'+str(obs_ct)) * xsapec('apec'+str(obs_ct))  # get model
    # define model parameters
    get_model_component('apec' + str(obs_ct)).kT = temperature
    get_model_component('apec' + str(obs_ct)).redshift = redshift  
    get_model_component('apec' + str(obs_ct)).Abundanc = abundance
    get_model_component('abs1').nH = nH_val  
    set_source(obs_count, src_model_dict[obsid]) #set model to source
    x = plot_model()
    print(x)
    return mdl

def Fitting(base_directory,dir,file_name,num_files,redshift,n_H,output_file, AGN):
    """
    Fit each region's spectra and create a text file containing the spectral
    fit information for each bin

    Args:
        base_directory: Path to main Directory
        dir: ObsID
        file_name: Root name of PI/PHA file
        num_files: Number of spatial bins
        redshift: Redshift of cluster
        n_H: Column density
        output_file: Text file containing each bin's spectral fit information

    Return:
        None
    """
    plot_dir = base_directory+'/FitPlots/'
    output_file = output_file.split('.')[0]
    os.chdir(base_directory)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    # Non Deprojected Fits
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature Temp_min Temp_max Abundance Ab_min Ab_max Norm Norm_min Norm_max ReducedChiSquare Flux Flux_min Flux_max\n")
    for i in range(num_files):
        print("Fitting model to spectrum number "+str(i+1))
        spectrum_files = []
        background_files = []
        model_list = []
        data_list = []
        
        for directory in dir:
            if num_files == 1:  # Are we fitting multiple files or not?
                spectrum_files.append(directory+'/repro/'+file_name+".pi")
                background_files.append(directory+'/repro/'+file_name+"_bkg.pi")
            else:
                spectrum_files.append(directory+'/repro/'+file_name+str(i)+".pi")
                background_files.append(directory+'/repro/'+file_name+str(i)+"_bkg.pi")
            data_list.append(load_pha(i,spectrum_files[i],use_errors=True)) #Read in)
            model_list.append(sherpa_model())                 
        # for each region calculate the model
        #     
    file_to_write.close()