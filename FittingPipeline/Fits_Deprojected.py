'''
65;6003;1c------------------------------------------------------
GOAL:
    Step through bins (spectra) and calculate the temperature value
    of each bin using XSPEC for DEPROJECTED quantites
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
from LSCalc import ls_calc
import deproject
from astropy import units as u
from sherpa.all import *
from sherpa.astro.all import *
from sherpa.astro.ui import *

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

def Fitting_Deprojected(base_directory,ObsIDs,file_name,num_files,redshift,n_H,Temp_guess,output_file, region_dir, reg_file_prefix, num_bins):
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
    plot_dir = base_directory+'/FitPlots/'
    output_file = output_file.split('.')[0]
    os.chdir(base_directory)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    # Set up deprojection
    # Get region values
    region_vals = []  # First get regions in kpc
    for region_ct in range(num_bins):
        with open(region_dir+reg_file_prefix+str(region_ct)+'.reg') as reg_:
            reg_data = reg_.readlines()[3].split(')')[0].split('(')[1]
            r_in_ = ls_calc(redshift,float(reg_data.split(',')[2].strip('"')))
            #r_in.append(r_in_)
            r_out_ = ls_calc(redshift,float(reg_data.split(',')[3].strip('"')))
            if r_in_ not in region_vals:
                region_vals.append(r_in_)
            if r_out_ not in region_vals:
                region_vals.append(r_out_)
    # Create deprojection instance
    dep = deproject.Deproject(radii=[float(x) for x in region_vals]* u.arcsec)
    # load associated datasets
    for ann in range(len(region_vals)-1):
        for obsid in ObsIDs:
            # Take pi file from obsid/repro
            dep.load_pha(obsid+'/repro/'+file_name+'%s.pi' % (str(ann)), annulus=ann)
    # Set up fit parameters
    dep.set_source('xsphabs*xsapec')
    dep.ignore(None, 0.5)
    dep.ignore(7.0, None)
    dep.freeze("xsphabs.nh")
    dep.set_par('xsapec.redshift', redshift)
    dep.set_par('xsphabs.nh', n_H)
    dep.set_par('xsapec.Abundanc', 0.4)
    dep.thaw('xsapec.Abundanc')
    set_method("levmar")
    set_stat("chi2xspecvar")
    dep.subtract()  # Subtract associated background. Read in automatically earlier
    onion = dep.fit()
    file_to_write = open(output_file+"_deproj.txt",'w+')
    file_to_write.write("BinNumber r_in r_out temp temp_min temp_max norm norm_min norm_max dens dens_min dens_max \n")
    #file_to_write.write('BinNumber r_in r_out temp norm dens \n')
    #try:
    onion_errs = dep.conf()
    for annulus in range(num_bins):
        #file_to_write.write('%i %.2E %.2E %.2E %.2E %.2E \n'%(onion['annulus'][annulus], onion['rlo_ang'][annulus], onion['rhi_ang'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.norm'][annulus], onion['density'][annulus]))

        file_to_write.write('%i %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E\n'%(onion['annulus'][annulus], onion['rlo_ang'][annulus], onion['rhi_ang'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.kT'][annulus]+onion_errs['xsapec.kT_lo'][annulus], onion['xsapec.kT'][annulus]+onion_errs['xsapec.kT_hi'][annulus], onion['xsapec.norm'][annulus], onion['xsapec.norm'][annulus]+onion_errs['xsapec.norm_lo'][annulus], onion['xsapec.norm'][annulus]+onion_errs['xsapec.norm_hi'][annulus], onion['density'][annulus], onion['density'][annulus]+onion_errs['density_lo'][annulus], onion['density'][annulus]+onion_errs['density_hi'][annulus]))
    
    """except:
        print('No error estimates due to a bad fit')
        for annulus in range(num_bins):
            #file_to_write.write('%i %.2E %.2E %.2E %.2E %.2E \n'%(onion['annulus'][annulus], onion['rlo_ang'][annulus], onion['rhi_ang'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.norm'][annulus], onion['density'][annulus]))
    
            file_to_write.write('%i %.2E %.2E %.2E %.2E %.2E \n'%(onion['annulus'][annulus], onion['rlo_ang'][annulus], onion['rhi_ang'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.kT'][annulus], onion['xsapec.norm'][annulus], onion['xsapec.norm'][annulus], onion['xsapec.norm'][annulus], onion['density'][annulus], onion['density'][annulus], onion['density'][annulus]))"""
    
    
        #file_to_write.write('%i %.2E %.2E %.2E %.2E %.2E \n'%(row[0], row[1], row[2], row[-3], row[-2], row[-1]))
    file_to_write.close()
