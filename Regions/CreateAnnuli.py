"""
This code will create annuli of the same SNR given an X-ray image

The X-ray image should **NOT** be exposure corrected!

"""
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from ciao_contrib.runtool import *


def calculate_annuli(image_path, center, targetSNR, num_annuli, output_dir, pix2arcsec=0.492):
    """
    Calculate the size of annuli containing the same SNR.

    Args:
        image_path: Path to image
        center: Center of the annuli in image coordiantes
        targetSNR: Target SNR
        num_annuli: Number of annuli
        output_dir: Output directory
        pix2arcsec: Pixel size in arcsec (default 0.492 for Chandra)

    Returns:
        inner_radii: List of inner radii
        outer_radii: List of outer radii

    """
    image = fits.getdata(image_path).T  # Get image data
    radius_inner = 0  # Inner radius in pixels
    radius_outer = 1  # Outer radius in pixels
    inner_radii = []  # List of inner radii
    outer_radii = []  # List of outer radii
    snr_list = []  # List of SNR
    cen_x, cen_y = center  # Obtain x and y coordinates

    annulus = 0  # Annulus number
    while annulus <= num_annuli:  # While the annulus number is less than the total number of desired annuli
        sum_flux = np.sum(
            np.sum(image[cen_x + radius_inner:cen_x + radius_outer, cen_y + radius_inner:cen_y + radius_outer], axis=0),
            axis=0)  # Calculate SNR
        SNR = np.sqrt(sum_flux)
        if SNR > targetSNR:
            annulus += 1  # Update annulus number
            inner_radii.append(radius_inner)  # Update inner radii list
            outer_radii.append(radius_outer)  # Update outer annuli list
            radius_inner = radius_outer  # Update inner radius
            snr_list.append(SNR)  # Update SNR list
        radius_outer += 1  # Update outer annulus
    x_values = [val * pix2arcsec for val in outer_radii]
    plt.plot(x_values, snr_list, 'o')
    plt.xlabel('distance from center (arcsec)', fontsize=16, fontweight='bold')
    plt.ylabel('signal-to-noise', fontsize=16, fontweight='bold')
    plt.ylim(0, 1.25 * np.max(snr_list))
    plt.xlim(x_values[0] - 5, x_values[-1] + 5)
    plt.hlines(targetSNR, xmin=x_values[0] - 10, xmax=x_values[-1] + 10, label='SNR Target', linestyle='--', color='g')
    plt.legend(loc='best')
    plt.savefig('%s/SNRvsDist.png' % output_dir)
    return inner_radii, outer_radii


def create_annuli_file(image_path, center, inner_radii, outer_radii, output_dir, pix2arcsec=0.492):
    """
    Create files containing the DS9 annuli. Each annulus will have its own file.

    Args:
        image_path: Path to image
        center: Center of the annuli in image coordiantes
        inner_radii: List of inner radii
        outer_radii: List of outer radii
        output_dir: Output directory
        pix2arcsec: Pixel size in arcsec (default 0.492 for Chandra)

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
    ra, dec = dmcoords.ra, dmcoords.dec
    annulus_ct = 0
    for inner_radius, outer_radius in zip(inner_radii, outer_radii):
        with open('%s/annulus_%i.reg' % (output_dir, annulus_ct), 'w+') as new_reg:
            new_reg.write("# Region file format: DS9 version 4.1 \n")
            new_reg.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
            new_reg.write("physical \n")
            #new_reg.write('annulus(%s,%s,%.3f",%.3f")' % (ra, dec, inner_radius * pix2arcsec, outer_radius * pix2arcsec))
            new_reg.write('annulus(%s,%s,%.3f,%.3f)' % (x, y, inner_radius * pix2arcsec, outer_radius * pix2arcsec))
        annulus_ct += 1


if __name__ == '__main__':
    image_path = '/home/crhea/Documents/PKS0745/PKS0745/0.5-8.0_thresh.img'
    center = (2466, 2194)
    output_dir = '/home/crhea/Documents/PKS0745/regions'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    inner_radii, outer_radii = calculate_annuli(image_path, center, 35, 10, output_dir)
    create_annuli_file(image_path, center, inner_radii, outer_radii, output_dir)
