base_dir = /home/crhea/Documents/PKS0745/Broad
#----------------------------WVT-----------------------------#
image_fits = broad_thresh.img
exposure_map = none
stn_target = 10
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_source_stn10
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/crhea/Documents/PKS0745/Broad
Name =
ObsIDs = 2427,12881
source_file = broad_thresh
output_dir = binned_source_10/
Temp_data = Temp_bin_source_stn10.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = True
bin_spec = True
num_bins = 486
fit_spec = True
plot = True
Colormap = inferno
