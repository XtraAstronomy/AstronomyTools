'''
Goal:
Create binned spectra from Chandra data given the WVT map of the pixels

OUTPUTS:
    - A combined spectra for each bin as designated by the WVT.

    - This is to be used for spectral fitting (we'll that's why I made this program)

    - File put in /PathToChandraData/OBSID/repro/binned

Additional Notes:
    As mentioned, the program was designed to generate combinned-binned-spectra
    so that I could generate temperature maps...
    The program can easily be canabilized for other uses or specifications
'''
import numpy as np
import os
import time
from ciao_contrib.runtool import *
from crates_contrib.utils import *
import shutil
from joblib import Parallel, delayed
from tqdm import tqdm
#-------------------------------------------------#
#-------------------------------------------------#
# Bin Class Information
# Everything in here should be self-explanatory... if not let me know and I
# will most happily comment it! :)
class Bin:
    def __init__(self, number):
        self.bin_number = number
        self.pixels = []
        self.total_pixels = 0
    def add_pixel(self, Pixel):
        self.pixels.append(Pixel)
        self.total_pixels += 1
#-------------------------------------------------#
#-------------------------------------------------#
# Pixel Class Information
# Ditto à propos the documentation for this class
class Pixel:
    def __init__(self, number, pix_x, pix_y, width, height):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
        self.width = width
        self.height = height
    def update(self, pix_x, pix_y, width, height):
        self.pix_x = pix_x
        self.pix_y = pix_y
        self.width = width
        self.height = height
#-------------------------------------------------#
#-------------------------------------------------#
# Get necessary filenames such as evt2,asol1,bpix1, and msk1
#   parameters:
#
def get_filenames(dir):
    filenames = dict()
    for file in os.listdir(dir):
        if file.endswith("evt2.fits"):
            filenames['evt2'] = file
        if file.endswith("bpix1.fits"):
            filenames['bpix1'] = file
    return filenames
#---------------------------------------------------#
#---------------------------------------------------#
# specextract_run
# apply specextract functio from ciao
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_convert = fits file in string format
#       outfile_from_convert = pha outroot in string format
#       output_dir = directory for output file
def specextract_run(obsid,filenames,file_to_convert,outfile_from_convert,output_dir,bin_number):
    """
    Execute specectract command with designated parameters. Grouping is set to 1 count per bin.

    Args:
        obsid (str): ObsID
        filenames (str): List of relavent files for given ObsID -- contains badpixel file, evt1, evt2, ...
        file_to_convert (str): Initial Level 2 file
        outfile_from_convert (str): Name of extracted spectrum
        output_dir (str): Directory for extracted spectrum
        bin_num (int): Relative number of the bin wrt WVT numbering system

    Returns:
        Creates extracted spectrum from WVT bin region
    """
    specextract.punlearn()
    specextract.infile = file_to_convert+"[sky=region("+str(output_dir)+str(bin_number)+"temp.reg)]" #1 to make sure we have enough space. serious overkill!
    specextract.outroot = outfile_from_convert
    specextract.bkgfile = obsid+'_blank.evt[sky=region('+str(output_dir)+str(bin_number)+'temp.reg)]'
    specextract.bkgresp = False #Necessary if using blank sky
    specextract.badpixfile = filenames['bpix1']
    specextract.grouptype = 'NUM_CTS'
    specextract.binspec = 1
    specextract.clobber = True
    specextract.energy_wmap = '500:14000'
    #print("     Running Specextract...")
    specextract()
    return True
#---------------------------------------------------#
#---------------------------------------------------#
#create temporary region file with regions in bin
#   parameters:
#       output_dir = output directory
#       reigons = list of regions to add
def create_reg(output_dir,regions):
    """
    Create temporary region file in ds9 format

    Args:
        output_dir (str): Output directory
        reigons (str): List of regions to add
    """
    #Create temporary reg file for split
    print(str(output_dir)+str(bin_number)+"temp.reg")
    with open(str(output_dir)+str(bin_number)+"temp.reg","w+") as file:
        file.write("# Region file format: DS9 version 4.1 \n")
        file.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
        file.write("physical \n")
        count = 0
        for region in regions:
            if count == 0:
                file.write(region+" \n")
            else:
                file.write('+'+region+" \n")
            count += 1
#---------------------------------------------------#
#---------------------------------------------------#
#create event file with regions in bin
#   parameters:
#
def create_evt(file_to_split,bin_number,output_dir):
    """
    Create temporary event file

    Args:
        file_to_split (str): Name of file to be split
        bin_number (int): Relative number of the bin wrt WVT numbering system
        output_dir (str): Path to output directory

    Returns:
        A copy of the temporary event file
    """
    dmcopy.punlearn()
    dmcopy.infile = file_to_split+'.fits[sky=region('+output_dir+str(bin_number)+'temp.reg)]'
    dmcopy.outfile = output_dir+'bin_'+str(bin_number)+'.fits'
    dmcopy.clobber = True
    dmcopy()
#---------------------------------------------------#
#---------------------------------------------------#
# Combine individuals pixels into larger boxes to
# reduce spexectract time.
def create_reg_comb(pix_in_bin,file_to_split,bin_number,output_dir):
    """
    Function to concatenate pixels in a bin together in order to reduce calculation time.

    Args:
        pix_in_bin (int): List of pixels in bin
        file_to_split (str): Name of event file
        bin_number (int): Relative number of the bin wrt WVT numbering system
        output_dir (str): Path to output directory

    Returns:
        Temporary ds9 file with concatenated pixels in bin that will be used to extract spectrum rapidly

    """
    points = {} #(x,[y]) pairs
    for pixel in pix_in_bin:
        if pixel.pix_x not in points.keys():
            points[pixel.pix_x] = [pixel.pix_y]
        else:
            points[pixel.pix_x].append(pixel.pix_y)
    pixel_x_list = [pixel.pix_x for pixel in pix_in_bin]
    pixel_y_list = [pixel.pix_y for pixel in pix_in_bin]
    x_min = np.min(pixel_x_list)
    y_min = np.min(pixel_y_list)
    x_max = np.max(pixel_x_list)
    y_max = np.max(pixel_y_list)
    pixels_used = {}

    for x in range(int(x_min),int(x_max)+1):
        if x not in pixels_used.keys():
            pixels_used[x] = []
    with open(output_dir+str(bin_number)+'temp.reg','w+') as file:#will create a temporary file used to create spectra
        file.write("# Region file format: CIAO version 1.0 \n")
        file_phys = open(output_dir+str(bin_number)+'temp_phys.reg','w+')
        file_phys.write("# Region file format: CIAO version 1.0 \n")
        #Step through each value to build up rectangle!
        in_count = 0
        for y in range(int(y_min),int(y_max)+1):
            for x in range(int(x_min),int(x_max)+1):
                xwi = 1 #Let us try and add a pixel to our box
                ywi = 1
                if y in points[x]: #only check point if it is in our bin
                    new_box = False #will we create a new region box?
                    continue_extension = True #we want to initialize the while loop for the box
                    if y in pixels_used[x]: #if we already assigned the pixel then dont contine
                        continue_extension = False
                        new_box = False
                    if y not in pixels_used[x]: #if not assigned yet we should assign
                        pixels_used[x].append(y)
                        new_box = True #that means we need a new box!
                    in_count += 1
                    while continue_extension == True: #Lets see how large we can make this box
                        if x+xwi > x_max: #Are we within bounds of bin?
                            break
                        if y+ywi > y_max: #ditto wrt y...
                            break
                        yextend = True
                        for xi in range(x,x+xwi+1): #can we add a y for all x-values currently used
                            if y+ywi in pixels_used[xi] or y+ywi not in points[xi]:
                                yextend = False
                                break
                        xextend = True
                        for yi in range(y,y+ywi+1): #ditto for x
                            if yi not in points[x+xwi] or yi in pixels_used[x+xwi]:
                                xextend = False
                                break
                        xandyextend = False #ok can we add the corner xandy extension point?
                        if x+xwi <= x_max: #lets make sure there is no conflicts
                            xandyextend = (y+ywi not in pixels_used[x+xwi] and y+ywi<=y_max)
                        if x+xwi > x_max:
                            pass
                        if (yextend and xextend and xandyextend): #ok were all good to extend our box!
                            for xi in range(x,x+xwi+1): #gotta get all those new x-values
                                if xi<=x_max:
                                    for yi in range(y,y+ywi+1): #and dont forget those y-values :)
                                        if yi not in pixels_used[xi]:
                                            pixels_used[xi].append(yi) #if it hasn't already been noted lets add it
                            xwi += 1 #ok now to make our box potentially bigger....
                            ywi += 1
                        if not (yextend and xextend and xandyextend):
                            continue_extension = False #fail!!!!
                    if new_box == True:#now we make a new box!
                        x_coord1,y_coord1,ra_center,dec_center = coord_trans((x+.5)+xwi*.5,(y+.5)+ywi*.5,file_to_split)
                        file.write("box(%.16f,%.16f,%i,%i) \n"%(x_coord1,y_coord1,xwi,ywi)) #Have to subtract 0.25 from x to make it align -- non permanent
                        file_phys.write("box(%s,%s,%i,%i) \n"%(ra_center,dec_center,xwi,ywi))
                else:
                    #value not in list of pixels
                    pass
        file_phys.close()
    shutil.copy(output_dir+str(bin_number)+'temp_phys.reg',output_dir+str(bin_number)+'.reg')
    return None
#---------------------------------------------------#
#---------------------------------------------------#
# Split up fits files into pi/pha files
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_split = fits file in string format
#       output_file = directory for output
#       output_split = pha outroot in string format
#       x_center = pixel x-position center in physical (sky) coordinates
#       y_center = pixel y-position center in physical (sky) coordinates
#       width = width of pixel in sky coordinates
#       height = height of pixel in sky coordinates
#       pix_in_bin_num = pixel number relative to bin (0-max(bin.pixels))
#       bin_number = Relative number of the bin wrt WVT numbering system
def split_fits(obsid,filenames,file_to_split,output_file,output_dir,pix_in_bin_num,bin_number):
    """
    Create spectra from initial bin regions

    Args:
        obsid (str): Current ObsID
        filenames (str): List of files necessary for specextract
        file_to_split (str): Fits file in string format
        output_file (str): Directory for output
        output_dir (str): Pha outroot in string format
        pix_in_bin_num (int): pixel number relative to bin (0-max(bin.pixels))
        bin_number (int): Relative number of the bin wrt WVT numbering system
    """
    output =  output_file+"_"+str(bin_number)+"_"
    regions = []
    file_to_convert = file_to_split+'.fits'
    create_reg_comb(pix_in_bin_num,file_to_split,bin_number,output_dir)
    create_evt(file_to_split,bin_number,output_dir)
    specextract_run(obsid,filenames,file_to_convert,output,output_dir,bin_number)
    return None
#---------------------------------------------------#
#---------------------------------------------------#
#Change from image or logical coordinate system to physical or sky coordinate system
#   parameters:
#       pixel_x  = pixels x coordinate in image coordinates
#       pixel_y  = pixels y coordinate in image coordinates
#       file_to_split = file which we will split later (needed to get new coordinates)
def coord_trans(pixel_x,pixel_y,file_to_split):
    """
    Translate image (logical) coordinates into physical (sky) coordinates

    Args:
        pixel_x (float): x coordinate in image coordinates
        pixel_y (float): y coordinate in image coordinates
        file_to_split (str): Name of event file

    Returns:
        x_center (float): x coordinate in physical coordinates
        y_center (float): y coordinate in physical coordinates
        ra (float): RA in degrees
        dec (float): DEC in degrees
    """

    tr = SimpleCoordTransform(file_to_split+'.img')
    x_center,y_center = tr.convert("image", "physical", pixel_x, pixel_y)
    dmcoords.punlearn()
    dmcoords.infile = file_to_split+'.img'  # OBSID+'_broad_thresh.img'
    dmcoords.option = 'logical'
    dmcoords.logicalx = pixel_x
    dmcoords.logicaly = pixel_y
    dmcoords()
    ra = dmcoords.ra
    dec = dmcoords.dec
    return x_center,y_center,ra,dec
#---------------------------------------------------#
def source_fits(filenames, source_file, obsid):
    '''
    Create fits and image file for primary source region in reprocessed Directory
    Also create blanksky file

    Args:
        filenames (str): dictionary containing evt2 file
        source_file (str): source region name without extension
        obsid (str): Chandra Observation ID

    Returns:
        Creates image file and blanksky file
    '''
    dmcopy.punlearn()
    dmcopy.infile = filenames['evt2']+'[sky=region('+source_file+'.reg)]'
    dmcopy.outfile = source_file+'.fits'
    dmcopy.clobber = True
    dmcopy()
    dmcopy.punlearn()
    dmcopy.infile = source_file+'.fits'
    dmcopy.outfile = source_file+'.img'
    dmcopy.option = 'IMAGE'
    dmcopy.clobber = True
    dmcopy()
    # Now for the blanksky background file
    blanksky.punlearn()
    blanksky.evtfile = filenames['evt2']
    blanksky.outfile = str(obsid)+'_blank.evt'
    blanksky.clobber = True
    blanksky()
    return None
#---------------------------------------------------#
def spec_loop(obsid,filenames,file_to_split,output_file,output_dir,directory_repro,bin_i):
    '''
    Parallelized loop for creating spectra for bins
    '''
    os.chdir(directory_repro)
    pix_in_bin = bin_i.pixels
    try:
    	split_fits(obsid,filenames,file_to_split,output_file,output_dir,pix_in_bin,bin_i.bin_number)
    except:
    	pass

    return None
#---------------------------------------------------#
# Main program to create bin pi/pha files
#   parameters:
#       base_directory = Directory containing Chandra data
#       filename = name of file to read in WVT bin data
#       dir = Directory for Chandra OBSID
#       file_to_split = Name of file to split in repro directory
#       output_dir = Output path for binned files
def create_spectra(base_directory,filename,OBSIDS,source_file,output_dir,wvt_output):
    """
    Wrapper function to create spectra for each bin in the WVT map

    Args:
        base_directory (str): Directory containing Chandra data
        filename (str): Name of file to read in WVT bin data
        dir (str): Directory for Chandra OBSID
        file_to_split (str): Name of file to split in repro directory
        output_dir (str): Output path for binned files

    Returns:
        Individual spectrum files (.pi) for source and background in each WVT bin

    Note:
        This is parallelized to run on 4 cores in order to speed up calculation time
    """
    print("Starting to bin spectra...")
    for obsid in OBSIDS:
        print("We are on obsid %s"%obsid)
        print("#-------------------------------------------------------------------------------------#")
        directory = base_directory+'/'+obsid+'/repro'
        # Copy region file into current directory
        shutil.copy(base_directory+'/regions/'+source_file+'.reg', directory)
        # Create fits and img file for primary region within current directory
        #Just making sure the directory exists and is empty :)
        if not os.path.exists(directory+"/"+output_dir):
            os.makedirs(directory+"/"+output_dir)
        if len(os.listdir(directory+"/"+output_dir)) != 0:
            print("Cleaning output directory of  files")
            for item in os.listdir(directory+"/"+output_dir):
                os.remove(os.path.join(directory+"/"+output_dir, item))
        os.chdir(directory)
        output_file = directory+'/'+output_dir+'/'+source_file
        file_to_split = directory+'/'+source_file
        #snagging some file names for later...
        filenames = get_filenames(directory)
        print(' Running Blanksky Background...')
        source_fits(filenames, source_file, obsid)
        bins = []
        number_bins = -1
        pix_num = 0
        #Lets get the data from our WVT file
        with open(base_directory+'/'+filename+'.txt') as f:
            next(f)
            next(f) #skip first two lines
            for line in f:
                if not int(line.split(" ")[2]) in [bins[i].bin_number for i in range(len(bins))]:
                    # create new Bin instance if doesnt exist already
                    new_bin = Bin(int(line.split(" ")[2]))
                    bins.append(new_bin)
                    number_bins += 1
                # Create new pixel
                new_pix = Pixel(pix_num,float(line.split(" ")[0]),float(line.split(" ")[1]),0,0) #set width and height to zero
                pix_num += 1
                # Add pixel to current bin
                bins[number_bins].add_pixel(new_pix)
        # Execute parallel spectral creation
        Parallel(n_jobs=4,prefer="processes")(delayed(spec_loop)(obsid,filenames,file_to_split,output_file,output_dir,directory,bin_) for bin_ in tqdm(bins))
    return len(bins)
