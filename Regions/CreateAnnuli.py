"""
This code will create annuli of the same SNR given an X-ray image

"""
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from ciao_contrib.runtool import *


def calculate_annuli(image, center, targetSNR, num_annuli, output_dir, pix2arcsec=0.492):
    """

    Args:
        image:
        center:
        targetSNR:
        num_annuli:

    Returns:

    """

    radius_inner = 0  # Inner radius in pixels
    radius_outer = 1  # Outer radius in pixels
    inner_radii = []  # List of inner radii
    outer_radii = []  # List of outer radii
    snr_list = []  # List of SNR
    cen_x, cen_y = center  # Obtain x and y coordinates

    annulus = 0  # Annulus number
    while annulus <= num_annuli:  # While the annulus number is less than the total number of desired annuli
        sum_flux = np.sum(np.sum(image[cen_x+radius_inner:cen_x+radius_outer, cen_y+radius_inner:cen_y+radius_outer], axis=0), axis=0)  # Calculate SNR
        SNR = np.sqrt(sum_flux)
        if SNR > targetSNR:
            annulus += 1  # Update annulus number
            inner_radii.append(radius_inner)  # Update inner radii list
            outer_radii.append(radius_outer)  # Update outer annuli list
            radius_inner = radius_outer  # Update inner radius
            snr_list.append(SNR)  # Update SNR list
        radius_outer += 1  # Update outer annulus
    plt.plot([val * pix2arcsec for val in outer_radii], snr_list, 'o')
    plt.xlabel('distance from center (arcsec)', fontsize=16, fontweight='bold')
    plt.ylabel('signal-to-noise', fontsize=16, fontweight='bold')
    plt.savefig('%s/SNRvsDist.png'%output_dir)
    return inner_radii, outer_radii

def create_annuli_file(image_path, center, inner_radii, outer_radii, output_dir):
    """

    Args:
        image:
        inner_radii:
        outer_radii:

    Returns:

    """
    cen_x, cen_y = center  # Coordinates in image units
    # calculate physical units from image units
    dmcoords.punlearn()
    dmcoords.infile = image_path
    dmcoords.option = 'logical'
    dmcoords.logicalx = cen_x
    dmcoords.logicaly = cen_y
    dmcoords()
    x, y = dmcoords.x, dmcoords.y
    annulus_ct = 0
    for inner_radius,outer_radius in zip(inner_radii, outer_radii):
        with open('%s/annulus_%i.reg'%(output_dir, annulus_ct), 'w+') as new_reg:
            new_reg.write("# Region file format: DS9 version 4.1 \n")
            new_reg.write("physical \n")
            new_reg.write('annulus(%s,%s,%.3f,%.3f)' % (x,y,inner_radius*0.492,outer_radius*0.492))
        annulus_ct += 1


image_path = '/home/crhea/Documents/Perseus/XrayFits/img_src.img'
image = fits.getdata(image_path)
image_hdr = fits.getheader(image_path)


center = (215, 144)
output_dir = '/home/crhea/Documents/Perseus/XrayFits/regions'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
inner_radii, outer_radii = calculate_annuli(image.T, center, 35, 10, output_dir)
create_annuli_file(image_path, center, inner_radii, outer_radii, output_dir)