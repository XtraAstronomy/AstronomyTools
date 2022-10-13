sFor the mechanics of what is going on, please see the the documentation :)

This program will results in:
- A Weighted Voronoi Tessellation Map
- Spectra for each bin
- Temperature and Abundance Maps
- Hardness Maps


In order to run this program you need the following:
1. Reprocessed Chandra ObsIDs --> you also need to merge them (`merge\_obs`). 
2. Fits image of region of interest (`roi.reg`) created from merged Observations. You can make this with the following command with `ciao` loaded:
    `dmcopy 'merged_evt.fits[sky=region(roi.reg)]' img_src.img opt=IMAGE`

On a note -- the region files have to be square (so use box in ds9).
**Your region file must be in base_dir+'regions/'!!

In order to clean the data, you can use `DataCleaning/DataCleaningPipeline.py`. It will also create the merged observations.

To run the program please supply those items (and other relevant info) in an input file (see Inputs folder for examples)
and execute the following command:

`python Temperature_Maps.py Inputs/example.i`
