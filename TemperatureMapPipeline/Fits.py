'''
GOAL:
    Step through bins (spectra) and calculate the temperature value
    of each bin using XSPEC

OUTPUTS:
    A file which contains the bin number and associated
    temperature and reduced chi-squared
ADDITIONAL NOTES:
    Be sure to run heainit first
'''
#from astropy.io import fits
import os
import subprocess
from sherpa.optmethods import LevMar
from sherpa.stats import LeastSq
from sherpa.plot import DataPlot
from sherpa.astro.xspec import *
from sherpa.astro.all import *
from sherpa.astro.ui import *
#from sherpa.all import *
from multiprocessing import Process, JoinableQueue
from joblib import Parallel, delayed
from tqdm import tqdm
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
def obsid_set(src_model_dict,bkg_model_dict,obsid, bkg_spec,obs_count,redshift,nH_val,Temp_guess):
    """ Set the observations' model using an absorbed apec model. Each observation will be tied
    with the other to allow simultaneous fitting.

    We additionally apply a background model.

    It contains two (frozen) thermal emission pieces:
        1. unabsorbed apec model with kT=0.18 keV
        2. absorbed bremmstrahlung emission with kT=40 keV


    Args:
        src_model_dict (dict): Dictionary of source models -- Created in FitXSPEC()
        bkg_model_dict (dict): Dictionary of background models -- Created in FitXSPEC()
        obsid (str): Current Chandra ObsID
        bkg_spec (str): Background file corresponding to ObsID
        obs_ct (int): Relative ObsID number
        redshift (float): Cosmological Redshift of Object
        nH_val (float): Column Density in direction of object -- units of 10^{-22} cm^{-2}
        Temp_guess (float): Estimate of Temperature Value

    Returns:
        Updated src_model_dict and bkg_model_dict

    """
    load_pha(obs_count,obsid) #Read in
    # If first observation ID
    if obs_count == 1:
        src_model_dict[obsid] = xsphabs('abs'+str(obs_count)) * xsapec('apec'+str(obs_count)) #set model and name
        # Change src model component values
        get_model_component('apec' + str(obs_count)).kT = Temp_guess
        get_model_component('apec' + str(obs_count)).redshift = redshift  # need to tie all together
        get_model_component('apec' + str(obs_count)).Abundanc = 0.3
        thaw(get_model_component('apec' + str(obs_count)).Abundanc)
        get_model_component('abs1').nH = nH_val  # change preset value
        freeze(get_model_component('abs1'))
    else:
        src_model_dict[obsid] = get_model_component('abs1') * xsapec('apec' + str(obs_count))
        get_model_component('apec'+str(obs_count)).kT = get_model_component('apec1').kT #link to first kT
        get_model_component('apec' + str(obs_count)).redshift = redshift
        get_model_component('apec' + str(obs_count)).Abundanc = get_model_component('apec1').Abundanc  # link to first kT

    # Setup background model
    bkg_model_dict[obsid] = xsapec('bkgApec'+str(obs_count))+get_model_component('abs1')*xsbremss('brem'+str(obs_count))
    set_bkg(obs_count, unpack_pha(bkg_spec))
    set_source(obs_count, src_model_dict[obsid]) #set model to source
    set_bkg_model(obs_count,bkg_model_dict[obsid])
    #Change bkg model component values
    get_model_component('bkgApec' + str(obs_count)).kT = 0.18
    freeze(get_model_component('bkgApec'+str(obs_count)).kT)
    get_model_component('brem' + str(obs_count)).kT = 40.0
    freeze(get_model_component('brem' + str(obs_count)).kT)
    return None
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def FitXSPEC(spectrum_files,background_files,redshift,n_H,temp_guess,grouping,spec_count,plot_dir, errors=False):
    """
    Fit observations simultaneously using models defined in obsid_set().

    Args:
        spectrum_files (list): List of spectrum files (i.e. ['obs1.pi', 'obs2.pi'])
        background_files (list): List of background files (i.e. ['obs1_bkg.pi', 'obs2_bkg.pi'])
        redshift (float): Cosmological Redshift of Object
        n_H (float): Column Density in direction of object -- units of 10^{-22} cm^{-2}
        Temp_guess (float): Estimate of Temperature Value
        grouping (int): Number of counts per bin (i.e. 5)
        spec_count (int): Local number of spectra with reference to the total number of spectra being fitted
        plot_dir (str): Path to directory containing plots

    Kwargs:
        errors (bool): Boolean to calculate 1-sigma errors (default=False)

    Returns:
        Temperature: Fitted Temperature Value in keV
        Temp_min: Minimum Fitted Temperature Value at 1-sigma in keV
        Temp_max: Maximum Fitted Temperature Value at 1-sigma in keV
        Abundance: Fitted Abundance Value in Z_{solar}
        Ab_min: Minimum Fitted Abundance Value at 1-sigma in Z_{solar}
        Ab_max: Maximum Fitted Abundance Value at 1-sigma in Z_{solar}
        Norm: Fitted Normalization Value
        Norm_min: Minimum Fitted Normalization Value at 1-sigma
        norm_max: Maximum Fitted Normalization Value at 1-sigma
        reduced_chi_sq: Reduced Chi Square Value for the Fit
        Flux: Fitted Flux Value in ergs/s/cm^{-2}
        Flux_min: Minimum Flux Value in ergs/s/cm^{-2}
        Flux_max: Maximum Flux Value in ergs/s/cm^{-2}

    Note:
        If the fits are taking too long, set errors=False. The min and max value will just be copies of the fitted value.

        Flux errors are always calculated

    """
    #print('Everything is set')
    set_stat('chi2gehrels')
    set_method('levmar')
    hdu_number = 1  #Want evnts so hdu_number = 1
    src_model_dict = {}; bkg_model_dict = {}
    obs_count = 1
    for spec_pha in spectrum_files:  # Set up model for spectra and background
        obsid_set(src_model_dict, bkg_model_dict, spec_pha, background_files[int(obs_count-1)], obs_count, redshift, n_H, temp_guess)
        obs_count += 1
    for ob_num in range(obs_count-1):  # Group and set range
        group_counts(ob_num+1,grouping)
        notice_id(ob_num+1,0.5,8.0)
    #print('Fitting')
    fit()  # FIT!
    #Get important values
    Temperature = apec1.kT.val
    Abundance = apec1.Abundanc.val;
    #Calculate norm as average value
    Norm = 0; Norm_min = 0; Norm_max = 0
    for id_ in range(len(spectrum_files)):
        Norm += get_model_component('apec'+str(id_+1)).norm.val #add up values
    Norm = Norm/len(spectrum_files)


    if errors == False:
        Temp_min = Temperature  # +mins[0]
        Temp_max = Temperature  # +maxes[0]
        Ab_min = Abundance  # +mins[1];
        Ab_max = Abundance  # +maxes[1]
        Norm_min = Norm  # +Norm_min/len(spectrum_files)
        Norm_max = Norm  # +Norm_max/len(spectrum_files)
    elif errors == True:
        set_covar_opt("sigma",1)
        #print('Calculating Errors')
        covar(get_model_component('apec1').kT,get_model_component('apec1').Abundanc)
        #----------Calculate min/max values---------#
        mins = list(get_covar_results().parmins)
        maxes = list(get_covar_results().parmaxes)
        for val in range(len(mins)):
            if isFloat(mins[val]) == False:
                mins[val] = 0.0
            if isFloat(maxes[val]) == False:
                maxes[val] = 0.0
            else:
                pass
        Temp_min = Temperature  # +mins[0]
        Temp_max = Temperature  # +maxes[0]
        Ab_min = Abundance  # +mins[1];
        Ab_max = Abundance  # +maxes[1]

        Norm = 0; Norm_min = 0; Norm_max = 0
        for id_ in range(len(spectrum_files)):
            Norm += get_model_component('apec'+str(id_+1)).norm.val  # add up values
            #get errors
            covar(get_model_component('apec'+str(id_+1)).norm)
            mins = list(get_covar_results().parmins)
            maxes = list(get_covar_results().parmaxes)
            for val in range(len(mins)):
                if isFloat(mins[val]) == False:
                    mins[val] = 0.0
                if isFloat(maxes[val]) == False:
                    maxes[val] = 0.0
                else:
                    pass
            Norm_min += mins[0]
            Norm_max += maxes[0]
        Norm_min = Norm + Norm_min/len(spectrum_files)
        Norm_max = Norm + Norm_max/len(spectrum_files)

    f = get_fit_results()
    reduced_chi_sq = f.rstat
    #---------Set up Flux Calculation----
    #print('Flux Calculations')
    #flux_calculation = sample_flux(get_model_component('apec1'), 0.01, 50.0, num=1000, confidence=68)[0]
    Flux = 0#flux_calculation[0]
    Flux_min = 0#flux_calculation[1]
    Flux_max = 0#flux_calculation[2]
    reset(get_model())
    reset(get_source())
    clean()

    return Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq, Flux, Flux_min, Flux_max

#---------------------------------------------------------#
def fit_loop(dir,bin_spec_dir,file_name,redshift,n_H,temp_guess,grouping,plot_dir,base_directory,i):
    """
    Loop for fitting spectra. This will simply call the FitXSPEC function and write out the results to a temporary file

    Args:
        bin_spec_dir: Path to extracted spectra for each bin within an ObsID
        dir (str): ObsID
        file_name (str): Root name of PI/PHA file
        redshift (float): Cosmological Redshift of Object
        n_H (float): Column Density in direction of object -- units of 10^{-22} cm^{-2}
        Temp_guess (float): Estimate of Temperature Value
        grouping (int): Number of counts per bin (i.e. 5)
        plot_dir (str): Path to plots
        base_directory (str): Path to main Directory
        i (int): Relative spectrum number

    Return:
        None
    """
    #print('Starting fit')
    os.chdir(base_directory)
    spectrum_files = []
    background_files = []
    fit_bool = True
    for directory in dir:  # Step through each ObsID
        try:  # Collect spectrum if it exists
            spectrum_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+".pi")
            background_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+"_bkg.pi")
        except:
            break
    try:
        Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max =\
            FitXSPEC(spectrum_files,background_files,redshift,n_H,temp_guess,grouping,i,plot_dir)
    except:
        Temperature = Temp_min = Temp_max = Abundance = Ab_min = Ab_max = Norm = Norm_min = Norm_max = reduced_chi_sq = Flux = Flux_min = Flux_max = 0
    with open('temp_'+str(i)+'.txt','w+') as out_temp:
        out_temp.write("%i %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E\n"% \
                (i,Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max))
    return None
#--------------------------------------------------------------------#
def concat_temp_data(num_spec, output_file):
    """
    Concatenate temperature information. Used if running on several processors
    Args:
        num_spec: Number of spectra (or bins)
        output_file: File name of final concatenated temperature data
    """
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature Temp_min Temp_max Abundance Ab_min Ab_max Norm Norm_min Norm_max ReducedChiSquare Flux Flux_min Flux_max\n")
    # Get the data for each temporary file and then delete
    for spec_i in range(num_spec):
        with open('temp_'+str(spec_i)+'.txt', 'r') as temp_file:
            for line in temp_file.readlines():
                file_to_write.write(line)
    file_to_write.close()
    return None
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
def Fitting(base_directory,dir,file_name,num_files,redshift,n_H,temp_guess,output_file,bin_spec_dir):
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
        temp_guess: Initial temperature guess
        output_file: Text file containing each bin's spectral fit information
        bin_spec_dir: Path to extracted spectra for each bin within an ObsID

    Return:
        None
    """
    energy_min = 0.5
    energy_max = 8.0
    grouping = 10
    plot_dir = base_directory+'/FitPlots/'
    output_file = output_file.split('.')[0]
    os.chdir(base_directory)
    # Make sure plotting directory exists
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    Parallel(n_jobs=4,prefer="processes")(delayed(fit_loop)\
            (dir,bin_spec_dir,file_name,redshift,n_H,temp_guess,grouping,plot_dir,
                base_directory,i) for i in tqdm(range(num_files)))
    concat_temp_data(num_files, output_file)
