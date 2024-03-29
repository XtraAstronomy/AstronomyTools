# AstronomyTools
A set of astronomy tools developed in Professor Julie Hlavacek-Larrondo's lab.


## How to use this repository
Each subfolder has an individual ReadMe file describing what each program does and how to use them. You will also find ipython notebook in `Data Cleaning`, `Fitting Pipeline`, and `TemperatureMapPipeline`.

### Astrometry
A handful of simple routines to do basic astrometry with Chandra data

### Data Cleaning
A robust data cleaning pipeline for low-count (or not) Chandra data. This will create clean evt2 files and merged images.

### Fitting Pipeline
This pipeline fits Chandra observations in a given set of regions to create thermodynamic profiles

### GGF
This pipeline calculates Gaussian Gradient Filter images from Chandra observations

### TemperatureMapPipeline
As the name suggests, this subfolder contains the code to create thermodynamic maps using a WVT algorithm

### WVT
Here we have our implementation of the WVT algorithm


## API Documentation
Please go to our readthedocs page for the API documentation:

https://astronomytools.readthedocs.io/en/latest/


## Pumpkin
If you are interested in our work on machine learning for X-ray spectral analysis, please see [Pumpkin](https://github.com/XtraAstronomy/Pumpkin).


## Citation
If you use any tools in this package, please use the following citation

```
Rhea C. L., Hlavacek-Larrondo J., Prasow-Émond M., 2020, X-tra Astronomy
Tools Software Package, DOI 10.17605/OSF.IO/VWCKS
```