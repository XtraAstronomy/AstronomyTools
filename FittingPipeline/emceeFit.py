"""
Use Bayesian inference to calculate the standard thermodyanmic parameters in annuli.

To use this code, please update the parameters found in the INPUTS section. Then run `python emceeFit.py`
"""


import os
import numpy as np
from sherpa.astro.ui import *
import emcee 
import corner
import pickle 
import matplotlib.pyplot as plt
from multiprocessing import Pool
from LSCalc import ls_calc, ds_calc


##### INPUTS #####
set_stat('chi2gehrels')
ObsIDs = ['2427', '12881']
num_annuli = 10
redshift = 0.1028
nH_value = 0.4028

##### CODE #####

# Calculate comoving distance for density calculation
dist = ds_calc(redshift)
Da = dist[0]*3.086e21 #conversion to cm from kpc
# Parameters for electron density
const = 1e7*np.sqrt(4*np.pi)
red_dep = Da*(1+redshift)
r_mid = []
temps = []
temps_min = [] 
temps_max = [] 
abundances = [] 
abundances_min = []
abundances_max = []
densities = []
densities_min = []
densities_max = []

# Make sure that the output files are there
if not os.path.exists("FitPlots/CornerPlots"):
    os.makedirs("FitPlots/CornerPlots")
if not os.path.exists("emceeSamples"):
        os.makedirs("emceeSamples")

def log_probability(params, num_ObsID):
    """
    Log probability functoin
    """
    # Check parameter bounds
    if (params[0] > 10 or params[0] < 0.1 or
        params[1] > 1.3 or params[1] < 0.1 or
        params[2] > 1 or params[2] < 1e-9):
        return -np.inf
    else:
        # Update model parameters
        plasma1.kT = params[0]
        plasma1.abundanc = params[1]
        plasma1.norm = params[2]
    # Calculate the Sherpa statistic (e.g., chi-square)
    stat_summed = np.sum([calc_stat(i) for i in range(num_ObsID)])
    # Calculate the log-likelihood
    log_likelihood = -0.5 * (stat_summed)
    # Include any priors (if needed)
    log_prior = 0  # Assuming flat priors for simplicity
    return log_likelihood + log_prior

outfile = open('BayesianFits.txt','w+')
outfile.write("annulus,temperature,temperatureMin,temperatureMax,abundance,abundanceMin,abundanceMax,norm,normMin,normMax,density,densityMin,densityMax\n")
for annulus in range(num_annuli): # for each annulus fit model
    # Load your data
    for obs_ct, obsID in enumerate(ObsIDs):
        load_pha(obs_ct, '%s/repro/annulus_%i.pi'%(obsID, annulus))
        subtract(obs_ct)  # subtract background
        if annulus == 0:
            set_source(obs_ct, xsphabs.abs1 * (xsapec.plasma1 + xszphabs.zphabs1 *(xspowerlaw.pow1 + xsgaussian.gauss1)))
            gauss1.LineE = 6.4
            pow1.PhoIndex = 1.8
        else:
            set_source(obs_ct, xsphabs.abs1 * xsapec.plasma1)
    # Define the background (if any)
    #ui.load_bkg('')

    # Define the model (e.g., a thermal plasma model)
    #set_source(0, xsphabs.abs1 * xsapec.plasma1)

    # Set initial parameter values
    abs1.nh = nH_value
    plasma1.redshift = redshift
    # Freeze parameters that should not vary during fitting
    freeze(abs1.nh)
    freeze(plasma1.redshift)
    # Thaw abundance
    thaw(plasma1.Abundanc)

    fit()
    # Initial guesses for parameters
    initial = [plasma1.kT.val, plasma1.Abundanc.val, plasma1.norm.val]#, 0.102]
    # Number of dimensions (parameters)
    ndim = len(initial)

    # Number of walkers
    nwalkers = 500 #  500

    # Initial positions of walkers
    pos = initial + np.random.randn(nwalkers, ndim) * [1.0, 0.1, 1e-5]
    # Set up the MCMC sampler
    # Number of steps to run the MCMC
    nsteps = 200  # 200
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=len(ObsIDs), pool=pool)
        # Run the MCMC sampler
        sampler.run_mcmc(pos, nsteps, progress=True)

    # Get the samples from the sampler
    samples = sampler.get_chain(discard=20, flat=True)

    # Plot the corner plot
    figure = corner.corner(samples, labels=["Temperature", "Abundance", "Norm"])

    plt.savefig("FitPlots/CornerPlots/annulus_%s"%annulus)
    plt.clf()
    # Get the parameter estimates
    kT_mcmc, abund_mcmc, norm_mcmc = map(lambda v: [v[1], v[2]-v[1], v[1]-v[0]], 
                                            zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    temps.append(kT_mcmc[0])
    temps_min.append(kT_mcmc[1])
    temps_max.append(kT_mcmc[2])
    abundances.append(abund_mcmc[0])
    abundances_min.append(abund_mcmc[1])
    abundances_max.append(abund_mcmc[2])
    
    pickle.dump(samples, open("emceeSamples/samples_%i"%annulus, "wb"))
    # Calculate Radii
    # Read in radius data
    r_in = 0
    r_out = 0
    with open('regions/annulus_%i.reg'%annulus) as reg_:
        reg_data = reg_.readlines()[3].split(')')[0].split('(')[1]
        r_in = ls_calc(redshift,float(reg_data.split(',')[2].strip('"')))
        r_out = ls_calc(redshift,float(reg_data.split(',')[3].strip('"')))
        r_mid.append((r_in+r_out)/2)
    # Calculate volume
    dist_out = r_out*3.086e21
    dist_in = r_in*3.086e21
    vol = (4/3)*np.pi*(dist_out**3-dist_in**3)
    # Calculate electron density
    density = [const*red_dep*np.sqrt((1.2*norm)/vol) for norm in norm_mcmc]
    densities.append(density[0])
    densities_min.append(density[1])
    densities_max.append(density[2])
    outfile.write(f"{annulus},{kT_mcmc[0]},{kT_mcmc[1]},{kT_mcmc[2]}, {abund_mcmc[0]},{abund_mcmc[1]},{abund_mcmc[2]},{norm_mcmc[0]},{norm_mcmc[1]},{norm_mcmc[2]},{density[0]},{density[1]},{density[2]}\n")
    #reset(get_model())
    #reset(get_source())
    #clean()
    
    
    
# TEMPERATURE
plt.errorbar(r_mid,temps,yerr=[temps_min,temps_max], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel('Temperature (keV)')
plt.savefig('Temperatures.png')
plt.clf()
# Abundance
plt.errorbar(r_mid,abundances,yerr=[abundances_min,abundances_max], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel(r'Abundance (Z$_{\odot}$)')
plt.savefig('Abundances.png')
plt.clf()
# Abundance
plt.errorbar(r_mid,densities,yerr=[densities_min,densities_max], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel(r'Density cm$^{-3}$')
plt.yscale("log")
plt.savefig('Densities.png')
plt.clf()



