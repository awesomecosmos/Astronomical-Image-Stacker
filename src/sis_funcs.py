# -*- coding: utf-8 -*-

# This file contains all necessary packages and functions required for the 
# Sidereal Image Stacking of the images. This file is required to run the file 
# sidereal_img_stacker.py.
# Import it as: 
# from sis_funcs import *

###############################################################################
#-------------------SECTION ONE: IMPORTING PACKAGES---------------------------#
###############################################################################

# path-type packages
import os
import glob
from pathlib import Path

# basic Python packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Astropy packages
import astropy
import astroalign as aa
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.utils.data import get_pkg_data_filename

# ccdproc packages
import ccdproc as ccdp
from ccdproc import Combiner
from ccdproc import ImageFileCollection
from ccdproc.utils.sample_directory import sample_directory_with_files

# Misc. packages
import seaborn as sns
from datetime import datetime

###############################################################################
#--------------------SECTION ONE: HELPER FUNCTIONS----------------------------#
###############################################################################

def registered_img_writer(img,registered_img,SIS_stack_path):
    """
    This function writes the original FITS header to a FITS file of the 
    registered image.
    
    Parameters
    ----------
    img : str
        Filename of image from which we want to extract header.
    
    registered_img : array
        Registered image as an array.
    
    SIS_stack_path : WindowsPath object
        Path to directory where sidereally-stacked image is to be saved.
    
    Returns
    -------
    registered_img_ccd : astropy.nddata.ccddata.CCDData
        CCDData object of registered image.
    """
    img_hdul = fits.open(img)
    img_hdr1 = img_hdul[0].header
                
    run_filename = img_hdr1['RUN'].strip(' ')
    target_name = img_hdr1['FIELD'].strip(' ')
    exptime = img_hdr1['EXPTIME']
    filter_colour = img_hdr1['COLOUR'].strip(' ')
    obs_set = img_hdr1['SET'].strip(' ')
    chip_num = img_hdr1['CHIP']
                
    filename_to_write = "aligned-{}-{}-{}-{}-{}-{}.fit".format(run_filename,target_name,
                                                               exptime,filter_colour,
                                                               obs_set,chip_num)
        
    registered_img_ccd = CCDData(registered_img,unit='adu')
    registered_img_ccd.header = img_hdr1
    registered_img_ccd.write(SIS_stack_path/filename_to_write,overwrite=True)
    
    return registered_img_ccd

def aligned_img_stats(aligned_image_ccd_lst,plots_path):
    median_count = []
    mean_count = []
    for a_file in aligned_image_ccd_lst:
        hdu = CCDData(a_file,unit='adu')
        image_data = hdu.data.astype(float) 
        image_hdr = hdu.header
        
        median_count.append(np.median(a_file))
        mean_count.append(np.mean(a_file))
    
    min_count_for_median = np.min(median_count)
    min_count_for_mean = np.min(mean_count)
    max_count_for_median = np.max(median_count)
    max_count_for_mean = np.max(mean_count)
    
    plt.figure()
    plt.plot(mean_count, label='mean',color="palevioletred")
    plt.axhline(y=min_count_for_mean,linestyle='-',linewidth=0.5,color='blue',label='min mean {:.2f}'.format(min_count_for_mean),alpha=1)
    plt.axhline(y=max_count_for_mean,linestyle='-',linewidth=0.5,color='blue',label='max mean {:.2f}'.format(max_count_for_mean),alpha=1)
    plt.xlabel('Image number')
    plt.ylabel('Count (ADU)')
    plt.title('Mean pixel value for aligned images')
    plt.legend()
    plt.grid()
    plt.savefig(plots_path/"aligned_stats_mean.jpg",dpi=900)
    plt.show()

    plt.figure()
    plt.plot(median_count, label='median',color="darkviolet")
    plt.axhline(y=min_count_for_median,linestyle='-',linewidth=0.5,color='red',label='min median {:.2f}'.format(min_count_for_median),alpha=1)
    plt.axhline(y=max_count_for_median,linestyle='-',linewidth=0.5,color='red',label='max median {:.2f}'.format(max_count_for_median),alpha=1)               
    plt.xlabel('Image number')
    plt.ylabel('Count (ADU)')
    plt.title('Median pixel value for aligned images')
    plt.legend()
    plt.grid()
    plt.savefig(plots_path/"aligned_stats_median.jpg",dpi=900)
    plt.show()