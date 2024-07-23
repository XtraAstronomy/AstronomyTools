import pandas as pd
import matplotlib.pyplot as plt
from LSCalc import ls_calc

redshift = 0.1028

data = pd.read_csv("BayesianFits.txt")
r_mid = []
for annulus in range(10):
    r_in = 0
    r_out = 0
    with open('regions/annulus_%i.reg'%annulus) as reg_:
        reg_data = reg_.readlines()[3].split(')')[0].split('(')[1]
        r_in = ls_calc(redshift,float(reg_data.split(',')[2].strip('"')))
        r_out = ls_calc(redshift,float(reg_data.split(',')[3].strip('"')))
        r_mid.append((r_in+r_out)/2)
# TEMPERATURE
plt.errorbar(r_mid,data['temperature'],yerr=[data['temperatureMin'],data['temperatureMax']], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel('Temperature (keV)')
plt.savefig('Temperatures.png')
plt.clf()
# Abundance
plt.errorbar(r_mid,data['abundance'],yerr=[data['abundanceMin'],data['abundanceMax']], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel(r'Abundance (Z$_{\odot}$)')
plt.savefig('Abundances.png')
plt.clf()
# Abundance
plt.errorbar(r_mid,data['density'],yerr=[data['densityMin'],data['densityMax']], fmt='o')
plt.xlabel('Radius (kpc)')
plt.ylabel(r'Density cm$^{-3}$')
plt.yscale("log")
plt.savefig('Densities.png')
plt.clf()