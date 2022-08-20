Create profile for galaxy cluster by fitting spectra

This program will results in:
- Spectra for each region
- Temperature, Abundance, Cooling Time, Pressure, and Density Plots


In order to run this program you need the following:
1. Reprocessed Chandra ObsIDs
  example: 4636/repro/

To run the program please supply those items (and other relevant info) in an input file (see Inputs folder for examples)
and execute the following command:

`python Temperature_Profiles.py Inputs/example.i`

Region files should be in the regions sub-directory. They should be labeled by name_count
where count corresponds to their annulus number. For example, if we have two regions, ann_1 and ann_2,
then ann_2 would be the outer annulus and ann_1 would be the inner annulus.

What does the program do?

- Create extracted region pi files for each ObsID in their reprocessed (repro) folders.
- Fit normal and deprojected spectra in each annulus
- Create text files and plots of thermodynamic variables (temperature, abundance, normalization, density, pressure, cooling time)

### FAQ

**Q:** What if I want to fit regions that are not ordered?

**A:** The easiest solution is to simply rename your regions to have the same prefix. For example, let's say you have two regions: `north.reg` and `south.reg`. In order to use the provided code without changing how inputs are read, you would need to rename them to be `reg_1.reg` and `reg_2.reg`. It doesn't matter which region corresponds to which number.


**Q:** Can I fit deprojected regions?

**A:** Yes you can! Simply set `deproj = True` in the input file. Keep in mind that this should *ONLY* be used if the regions are concentric annuli! Note that the background is subtracted in the case of deprojected fitting.

**Q:** What kind of model is fit?

**A:** We fit an absorbed thermal model to the spectrum: `phabs*apec`. However, you can easily change this in the fitting routines. Our background model has two components: `apec + phabs*brem`. This model was taken from McDonald et al. 2010 (need to add citation). Again, this can easily be changed :)
