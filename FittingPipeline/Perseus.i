#--------------------------SPECTRUM EXTRACTION-----------------------------#
reg_file_prefix = annulus_
num_files = 10
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/carterrhea/carterrhea/Perseus/Test
Name =
ObsIDs = 3209,4289
source_file = reg
output_dir = binned/
Temp_data = Temp_bin.txt
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
extract_spectrum = False
fit_spec = True
deproj = False
plot = True
